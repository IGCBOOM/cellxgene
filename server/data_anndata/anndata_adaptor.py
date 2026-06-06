import os
import tempfile
import time
import warnings
import importlib.metadata
from concurrent.futures import ThreadPoolExecutor
from threading import RLock
from uuid import uuid4

import anndata
import numpy as np
import pandas as pd
from packaging import version
from pandas.core.dtypes.dtypes import CategoricalDtype
from scipy import sparse

import server.common.compute.diffexp_generic as diffexp_generic
import server.common.compute.estimate_distribution as estimate_distribution
from server.common.colors import convert_anndata_category_colors_to_cxg_category_colors
from server.common.constants import Axis, MAX_LAYOUTS, XApproximateDistribution
from server.common.corpora import corpora_get_props_from_anndata
from server.common.errors import PrepareError, DatasetAccessError
from server.common.utils.type_conversion_utils import get_schema_type_hint_of_array
from server.data_common.data_adaptor import DataAdaptor
from server.common.fbs.matrix import encode_matrix_fbs
from server.common.compute.recluster import (
    OUTSIDE_RECLUSTER_CATEGORY,
    ReclusterJob,
    public_gene_filter_summary,
    recluster_from_expression_genes,
    recluster_from_obsm,
    resolve_gene_filter,
)

anndata_version = version.parse(str(importlib.metadata.version('anndata'))).release


def anndata_version_is_pre_070():
    major = anndata_version[0]
    minor = anndata_version[1] if len(anndata_version) > 1 else 0
    return major == 0 and minor < 7


class AnndataAdaptor(DataAdaptor):
    def __init__(self, data_locator, app_config=None, dataset_config=None):
        super().__init__(data_locator, app_config, dataset_config)
        self.data = None
        self.X_approximate_distribution = None
        self._recluster_jobs = {}
        self._recluster_results = {}
        self._recluster_lock = RLock()
        self._recluster_executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="cxg-recluster")
        self._load_data(data_locator)
        self._validate_and_initialize()

    def cleanup(self):
        self._recluster_executor.shutdown(wait=False, cancel_futures=True)

    @staticmethod
    def pre_load_validation(data_locator):
        if data_locator.islocal():
            # if data locator is local, apply file system conventions and other "cheap"
            # validation checks.  If a URI, defer until we actually fetch the data and
            # try to read it.  Many of these tests don't make sense for URIs (eg, extension-
            # based typing).
            if not data_locator.exists():
                raise DatasetAccessError("does not exist")
            if not data_locator.isfile():
                raise DatasetAccessError("is not a file")

    @staticmethod
    def file_size(data_locator):
        return data_locator.size() if data_locator.islocal() else 0

    @staticmethod
    def open(data_locator, app_config, dataset_config=None):
        return AnndataAdaptor(data_locator, app_config, dataset_config)

    def get_corpora_props(self):
        return corpora_get_props_from_anndata(self.data)

    def get_name(self):
        return "cellxgene anndata adaptor version"

    def get_library_versions(self):
        versions = dict(anndata=str(importlib.metadata.version("anndata")))
        for package_name in ("scanpy", "igraph", "leidenalg", "umap-learn"):
            try:
                versions[package_name] = str(importlib.metadata.version(package_name))
            except importlib.metadata.PackageNotFoundError:
                pass
        return versions

    @staticmethod
    def _create_unique_column_name(df, col_name_prefix):
        """given the columns of a dataframe, and a name prefix, return a column name which
        does not exist in the dataframe, AND which is prefixed by `prefix`

        The approach is to append a numeric suffix, starting at zero and increasing by
        one, until an unused name is found (eg, prefix_0, prefix_1, ...).
        """
        suffix = 0
        while f"{col_name_prefix}{suffix}" in df:
            suffix += 1
        return f"{col_name_prefix}{suffix}"

    def _alias_annotation_names(self):
        """
        The front-end relies on the existance of a unique, human-readable
        index for obs & var (eg, var is typically gene name, obs the cell name).
        The user can specify these via the --obs-names and --var-names config.
        If they are not specified, use the existing index to create them, giving
        the resulting column a unique name (eg, "name").

        In both cases, enforce that the result is unique, and communicate the
        index column name to the front-end via the obs_names and var_names config
        (which is incorporated into the schema).
        """
        self.original_obs_index = self.data.obs.index

        for ax_name, var_name in ((Axis.OBS, "obs"), (Axis.VAR, "var")):
            config_name = f"single_dataset__{var_name}_names"
            parameter_name = f"{var_name}_names"
            name = getattr(self.server_config, config_name)
            df_axis = getattr(self.data, str(ax_name))
            if name is None:
                # Default: create unique names from index
                if not df_axis.index.is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.index must be unique. "
                        "Please prepare data to contain unique index values, or specify an "
                        "alternative with --{ax_name}-name."
                    )
                name = self._create_unique_column_name(df_axis.columns, "name_")
                self.parameters[parameter_name] = name
                # reset index to simple range; alias name to point at the
                # previously specified index.
                df_axis.rename_axis(name, inplace=True)
                df_axis.reset_index(inplace=True)
            elif name in df_axis.columns:
                # User has specified alternative column for unique names, and it exists
                if not df_axis[name].is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.{name} must be unique. " "Please prepare data to contain unique values."
                    )
                df_axis.reset_index(drop=True, inplace=True)
                self.parameters[parameter_name] = name
            else:
                # user specified a non-existent column name
                raise KeyError(f"Annotation name {name}, specified in --{ax_name}-name does not exist.")

    def _create_schema(self):
        self.schema = {
            "dataframe": {
                "nObs": self.cell_count,
                "nVar": self.gene_count,
                **get_schema_type_hint_of_array(self.data.X),
            },
            "annotations": {
                "obs": {"index": self.parameters.get("obs_names"), "columns": []},
                "var": {"index": self.parameters.get("var_names"), "columns": []},
            },
            "layout": {"obs": []},
        }
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                ann_schema = {"name": ann, "writable": False}
                ann_schema.update(get_schema_type_hint_of_array(curr_axis[ann]))
                self.schema["annotations"][ax]["columns"].append(ann_schema)

        for layout in self.get_embedding_names():
            layout_schema = {"name": layout, "type": "float32", "dims": [f"{layout}_0", f"{layout}_1"]}
            self.schema["layout"]["obs"].append(layout_schema)

    def get_schema(self):
        return self.schema

    def _load_data(self, data_locator):
        # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
        # cost of significantly slower access to X data.
        try:
            # there is no guarantee data_locator indicates a local file.  The AnnData
            # API will only consume local file objects.  If we get a non-local object,
            # make a copy in tmp, and delete it after we load into memory.
            with data_locator.local_handle() as lh:
                # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
                # cost of significantly slower access to X data.
                backed = "r" if self.server_config.adaptor__anndata_adaptor__backed else None
                self.data = anndata.read_h5ad(lh, backed=backed)

        except ValueError:
            raise DatasetAccessError(
                "File must be in the .h5ad format. Please read "
                "https://github.com/theislab/scanpy_usage/blob/master/170505_seurat/info_h5ad.md to "
                "learn more about this format. You may be able to convert your file into this format "
                "using `cellxgene prepare`, please run `cellxgene prepare --help` for more "
                "information."
            )
        except MemoryError:
            raise DatasetAccessError("Out of memory - file is too large for available memory.")
        except Exception as e:
            import traceback
            error_msg = str(e)

            # IMPROVEMENT: Broadly catch ANY version incompatibility
            if "No read method registered" in error_msg and "IOSpec" in error_msg:
                message = (
                    "Error loading file: This H5AD file uses a newer internal format that "
                    "your version of 'anndata' cannot read.\n"
                    f"The specific error was: {error_msg}\n"
                    "Please upgrade anndata in your environment (pip install --upgrade anndata)."
                )
            else:
                message = (
                    "File not found or is inaccessible. File must be an .h5ad object. "
                    "Please check your input and try again."
                )

            if self.server_config.app__verbose:
                message += f"\n{traceback.format_exc()}"
            raise DatasetAccessError(message)

    def _validate_and_initialize(self):
        if anndata_version_is_pre_070():
            warnings.warn(
                "Use of anndata versions older than 0.7 will have serious issues. Please update to at "
                "least anndata 0.7 or later."
            )

        # var and obs column names must be unique
        if not self.data.obs.columns.is_unique or not self.data.var.columns.is_unique:
            raise KeyError("All annotation column names must be unique.")

        self._alias_annotation_names()
        self._validate_data_types()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self._create_schema()
        self.parameters.update({
            "recluster-enabled": True,
            "recluster-representations": self.get_recluster_representations(),
        })

        if self.dataset_config.X_approximate_distribution == "auto":
            """Lazy evaluate the heuristic if we are backed."""
            if not self.data.isbacked:
                self.X_approximate_distribution = estimate_distribution.estimate_approximate_distribution(self.data.X)
        else:
            self.X_approximate_distribution = self.dataset_config.X_approximate_distribution

        # heuristic
        n_values = self.data.shape[0] * self.data.shape[1]
        if (n_values > 1e8 and self.server_config.adaptor__anndata_adaptor__backed is True) or (n_values > 5e8):
            self.parameters.update({"diffexp-may-be-slow": True})


    def get_recluster_representations(self):
        """Return .obsm keys suitable for memory-conscious reclustering."""
        reps = []
        for key, value in self.data.obsm.items():
            if not isinstance(key, str):
                continue
            try:
                if len(value.shape) == 2 and value.shape[0] == self.data.n_obs and value.shape[1] >= 2:
                    reps.append(key)
            except Exception:
                continue
        if "X_pca" in reps:
            reps.remove("X_pca")
            reps.insert(0, "X_pca")
        return reps

    def get_recluster_gene_names(self):
        """Return frontend-visible gene names used by the gene-list UI.

        CELLxGENE may alias var.index into a unique annotation column during
        load. Use that configured display/index column so pasted gene symbols
        match what users see in the interface.
        """
        var_name_col = self.parameters.get("var_names")
        if var_name_col is not None and var_name_col in self.data.var:
            return self.data.var[var_name_col].astype(str).to_numpy()
        return np.asarray([str(x) for x in self.data.var_names], dtype=object)

    def _resolve_recluster_gene_filter(self, params):
        return resolve_gene_filter(self.get_recluster_gene_names(), params)

    def _validate_recluster_gene_filter_limits(self, n_obs, gene_filter):
        if gene_filter is None:
            return

        n_vars = int(gene_filter["included_gene_count"])
        if self.server_config.exceeds_limit("recluster_gene_count_max", n_vars):
            raise ValueError(
                f"Expression-gene reclustering would use {n_vars} genes, which exceeds "
                f"the configured recluster_gene_count_max limit"
            )

        expression_values = int(n_obs) * n_vars
        if self.server_config.exceeds_limit("recluster_expression_values_max", expression_values):
            raise ValueError(
                "Expression-gene reclustering would materialize too many expression values "
                f"({expression_values}). Use a smaller cell selection, fewer genes, "
                "or raise recluster_expression_values_max in the server config."
            )

    def _recluster_cleanup_locked(self):
        """Evict old completed jobs/results to bound memory use."""
        ttl = self.server_config.limits__recluster_result_ttl_seconds
        now = time.time()
        if ttl is not None:
            stale_result_ids = [
                result_id
                for result_id, result in self._recluster_results.items()
                if now - result.created_at > ttl
            ]
            for result_id in stale_result_ids:
                self._recluster_results.pop(result_id, None)

            stale_job_ids = [
                job_id
                for job_id, job in self._recluster_jobs.items()
                if job.status in ("complete", "error") and now - job.updated_at > ttl
            ]
            for job_id in stale_job_ids:
                self._recluster_jobs.pop(job_id, None)

        max_results = self.server_config.limits__recluster_results_per_session
        if max_results is not None and max_results >= 0:
            results_by_user = {}
            for result in self._recluster_results.values():
                results_by_user.setdefault(result.user_id, []).append(result)
            for user_id, results in results_by_user.items():
                results.sort(key=lambda r: r.created_at, reverse=True)
                for old in results[max_results:]:
                    self._recluster_results.pop(old.result_id, None)

    def _get_recluster_result(self, user_id, result_id):
        with self._recluster_lock:
            result = self._recluster_results.get(result_id)
            if result is None or result.user_id != user_id:
                raise KeyError(f"Unknown recluster result: {result_id}")
            return result

    def _update_recluster_job(self, job_id, *, status=None, stage=None, progress=None, result=None, error=None):
        with self._recluster_lock:
            job = self._recluster_jobs[job_id]
            job.update(status=status, stage=stage, progress=progress)
            if result is not None:
                job.result = result
            if error is not None:
                job.error = error

    def _run_recluster_job(self, job_id, user_id, obs_indices, params, result_id, gene_filter=None):
        def progress(stage, value):
            self._update_recluster_job(job_id, status="running", stage=stage, progress=value)

        try:
            progress("Preparing selected cells", 0.05)
            if gene_filter is None:
                result = recluster_from_obsm(
                    self.data,
                    obs_indices,
                    result_id=result_id,
                    user_id=user_id,
                    params=params,
                    progress=progress,
                )
            else:
                var_indices = gene_filter["var_indices"]
                gene_names = self.get_recluster_gene_names()[var_indices]
                result = recluster_from_expression_genes(
                    self.data,
                    obs_indices,
                    var_indices=var_indices,
                    var_names=gene_names,
                    gene_filter=gene_filter,
                    result_id=result_id,
                    user_id=user_id,
                    params=params,
                    progress=progress,
                )
            with self._recluster_lock:
                self._recluster_results[result.result_id] = result
                self._recluster_cleanup_locked()
            self._update_recluster_job(
                job_id,
                status="complete",
                stage="Complete",
                progress=1.0,
                result=result.schema_payload(),
            )
        except Exception as e:  # noqa: B902 - preserve the message for the polling UI
            self._update_recluster_job(
                job_id,
                status="error",
                stage="Error",
                progress=1.0,
                error=str(e),
            )

    def recluster_obs_start(self, user_id, args):
        args = args or {}
        try:
            obs_filter = args.get("filter", {}).get("obs", None)
            if obs_filter is None:
                raise ValueError("missing filter.obs")
            obs_mask = self._axis_filter_to_mask(Axis.OBS, obs_filter, self.data.n_obs)
        except (KeyError, IndexError, TypeError) as e:
            raise ValueError(f"Error parsing obs filter: {e}")

        obs_indices = np.flatnonzero(obs_mask).astype(np.int64)
        n_obs = int(obs_indices.size)
        if n_obs < self.server_config.limits__recluster_cellcount_min:
            raise ValueError(
                f"Select at least {self.server_config.limits__recluster_cellcount_min} cells to recluster"
            )
        if self.server_config.exceeds_limit("recluster_cellcount_max", n_obs):
            raise ValueError("Recluster request exceeds max cell count limit")

        params = dict(args.get("params", {}) or {})
        gene_filter = self._resolve_recluster_gene_filter(params)
        self._validate_recluster_gene_filter_limits(n_obs, gene_filter)

        if gene_filter is None:
            use_rep = params.get("use_rep", "X_pca")
            if use_rep not in self.data.obsm:
                available = ", ".join(self.get_recluster_representations())
                raise ValueError(f"Representation {use_rep!r} is not available for reclustering. Available: {available}")

        gene_filter_summary = public_gene_filter_summary(gene_filter)
        n_vars = None if gene_filter is None else int(gene_filter["included_gene_count"])

        with self._recluster_lock:
            self._recluster_cleanup_locked()
            max_jobs = self.server_config.limits__recluster_concurrent_jobs
            running_for_user = [
                job
                for job in self._recluster_jobs.values()
                if job.user_id == user_id and job.status in ("queued", "running")
            ]
            if max_jobs is not None and len(running_for_user) >= max_jobs:
                raise ValueError("A reclustering job is already running for this session")

            job_id = f"job_{uuid4().hex[:12]}"
            result_id = uuid4().hex[:12]
            job = ReclusterJob(
                job_id=job_id,
                user_id=user_id,
                n_obs=n_obs,
                n_vars=n_vars,
                gene_filter=gene_filter_summary,
            )
            self._recluster_jobs[job_id] = job

        self._recluster_executor.submit(
            self._run_recluster_job,
            job_id,
            user_id,
            obs_indices,
            params,
            result_id,
            gene_filter,
        )
        return job.to_dict()

    def recluster_obs_job_get(self, user_id, job_id):
        with self._recluster_lock:
            job = self._recluster_jobs.get(job_id)
            if job is None or job.user_id != user_id:
                raise KeyError(f"Unknown recluster job: {job_id}")
            return job.to_dict()

    def recluster_layout_to_fbs_matrix(self, user_id, result_id, layout_name):
        result = self._get_recluster_result(user_id, result_id)
        if layout_name != result.layout_name:
            raise KeyError(f"Unknown recluster layout: {layout_name}")

        full_embedding = np.full((self.data.n_obs, 2), np.nan, dtype=np.float32)
        full_embedding[result.obs_indices, :] = result.embedding[:, 0:2]
        normalized_layout = DataAdaptor.normalize_embedding(full_embedding)
        df = pd.DataFrame(normalized_layout, columns=[f"{layout_name}_0", f"{layout_name}_1"])
        return encode_matrix_fbs(df, col_idx=df.columns, row_idx=None)

    def recluster_annotation_to_fbs_matrix(self, user_id, result_id, annotation_name):
        result = self._get_recluster_result(user_id, result_id)
        if annotation_name != result.cluster_name:
            raise KeyError(f"Unknown recluster annotation: {annotation_name}")

        labels = np.full((self.data.n_obs,), OUTSIDE_RECLUSTER_CATEGORY, dtype=object)
        labels[result.obs_indices] = result.leiden
        categories = result.categories + [OUTSIDE_RECLUSTER_CATEGORY]
        series = pd.Series(pd.Categorical(labels, categories=categories), name=annotation_name)
        df = pd.DataFrame({annotation_name: series})
        return encode_matrix_fbs(df, col_idx=df.columns, row_idx=None)

    def recluster_export_h5ad(self, user_id, result_id):
        result = self._get_recluster_result(user_id, result_id)
        if self.data.isbacked:
            subset = self.data[result.obs_indices, :].to_memory()
        else:
            subset = self.data[result.obs_indices, :].copy()

        subset.obsm["X_cellxgene_recluster_umap"] = result.embedding
        subset.obs["cellxgene_recluster_leiden"] = pd.Categorical(result.leiden, categories=result.categories)
        subset.obs["cellxgene_original_obs_index"] = result.obs_indices
        if result.var_indices is not None:
            used_gene_mask = np.zeros((self.data.n_vars,), dtype=bool)
            used_gene_mask[result.var_indices] = True
            subset.var["cellxgene_recluster_gene_used"] = used_gene_mask
        subset.uns["cellxgene_recluster"] = {
            "result_id": result.result_id,
            "params": result.params,
            "n_obs": int(result.obs_indices.size),
            "n_vars": None if result.var_indices is None else int(result.var_indices.size),
        }

        fd, path = tempfile.mkstemp(prefix=f"cellxgene_{result_id}_", suffix=".h5ad")
        os.close(fd)
        subset.write_h5ad(path)
        return path, f"cellxgene_recluster_{result_id}.h5ad"

    def _is_valid_layout(self, arr):
        """return True if this layout data is a valid array for front-end presentation:
        * ndarray, dtype float/int/uint
        * with shape (n_obs, >= 2)
        * with all values finite or NaN (no +Inf or -Inf)
        """
        is_valid = type(arr) is np.ndarray and arr.dtype.kind in "fiu"
        is_valid = is_valid and arr.shape[0] == self.data.n_obs and arr.shape[1] >= 2
        is_valid = is_valid and not np.any(np.isinf(arr)) and not np.all(np.isnan(arr))
        return is_valid

    def _validate_data_types(self):
        # The backed API does not support interrogation of the underlying sparsity or sparse matrix type
        # Fake it by asking for a small subarray and testing it.   NOTE: if the user has ignored our
        # anndata <= 0.7 warning, opted for the --backed option, and specified a large, sparse dataset,
        # this "small" indexing request will load the entire X array. This is due to a bug in anndata<=0.7
        # which will load the entire X matrix to fullfill any slicing request if X is sparse.  See
        # user warning in _load_data().
        X0 = self.data.X[0, 0:1]
        if sparse.isspmatrix(X0) and not sparse.isspmatrix_csc(X0):
            warnings.warn(
                "Anndata data matrix is sparse, but not a CSC (columnar) matrix.  "
                "Performance may be improved by using CSC."
            )
        if self.data.X.dtype > np.dtype(np.float32):
            warnings.warn(
                f"Anndata data matrix is in {self.data.X.dtype} format not float32. " f"Precision may be truncated."
            )
        if self.data.X.dtype < np.float32:
            if self.data.isbacked:
                raise DatasetAccessError(
                    f"Data matrix in {self.data.X.dtype} format is not supported in backed mode."
                    " Please reload without --backed, or convert matrix to float32"
                )
            warnings.warn(
                f"Anndata data matrix is in unsupported {self.data.X.dtype} format -- will be cast to float32"
            )
            self.data.X = self.data.X.astype(np.float32)
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                datatype = curr_axis[ann].dtype
                downcast_map = {
                    "int64": "int32",
                    "uint32": "int32",
                    "uint64": "int32",
                    "float64": "float32",
                }
                if datatype in downcast_map:
                    warnings.warn(
                        f"Anndata annotation {ax}:{ann} is in unsupported format: {datatype}. "
                        f"Data will be downcast to {downcast_map[datatype]}."
                    )
                if isinstance(datatype, CategoricalDtype):
                    category_num = len(curr_axis[ann].dtype.categories)
                    if category_num > 500 and category_num > self.dataset_config.presentation__max_categories:
                        warnings.warn(
                            f"{str(ax).title()} annotation '{ann}' has {category_num} categories, this may be "
                            f"cumbersome or slow to display. We recommend setting the "
                            f"--max-category-items option to 500, this will hide categorical "
                            f"annotations with more than 500 categories in the UI"
                        )

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        if axis == Axis.OBS:
            if labels is not None and not labels.empty:
                df = self.data.obs.join(labels, self.parameters.get("obs_names"))
            else:
                df = self.data.obs
        else:
            df = self.data.var

        if fields is not None and len(fields) > 0:
            df = df[fields]
        return encode_matrix_fbs(df, col_idx=df.columns)

    def get_embedding_names(self):
        """
        Return pre-computed embeddings.

        function:
            a) generate list of default layouts
            b) validate layouts are legal.  remove/warn on any that are not
            c) cap total list of layouts at global const MAX_LAYOUTS
        """
        # load default layouts from the data.
        layouts = self.dataset_config.embeddings__names

        if layouts is None or len(layouts) == 0:
            layouts = [key[2:] for key in list(self.data.obsm.keys()) if type(key) is str and key.startswith("X_")]

        # remove invalid layouts
        valid_layouts = []
        obsm_keys = list(self.data.obsm.keys())
        for layout in layouts:
            layout_name = f"X_{layout}"
            if layout_name not in obsm_keys:
                warnings.warn(f"Ignoring unknown layout name: {layout}.")
            elif not self._is_valid_layout(self.data.obsm[layout_name]):
                warnings.warn(f"Ignoring layout due to malformed shape or data type: {layout}")
            else:
                valid_layouts.append(layout)

        if len(valid_layouts) == 0:
            raise PrepareError("No valid layout data.")

        # cap layouts to MAX_LAYOUTS
        return valid_layouts[0:MAX_LAYOUTS]

    def get_embedding_array(self, ename, dims=2):
        full_embedding = self.data.obsm[f"X_{ename}"]
        return full_embedding[:, 0:dims]

    def compute_diffexp_ttest(self, maskA, maskB, top_n=None, lfc_cutoff=None):
        if top_n is None:
            top_n = self.dataset_config.diffexp__top_n
        if lfc_cutoff is None:
            lfc_cutoff = self.dataset_config.diffexp__lfc_cutoff
        return diffexp_generic.diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff)

    def get_colors(self):
        return convert_anndata_category_colors_to_cxg_category_colors(self.data)

    def get_X_array(self, obs_mask=None, var_mask=None):
        # H5Py does not support boolean indexing (masks), so convert to integer indexing
        # when backed (ie, when AnnData is using H5Py indexing)
        if obs_mask is None:
            obs_mask = slice(None)
        elif self.data.isbacked and obs_mask.dtype == bool:
            obs_mask = obs_mask.nonzero()[0]
        if var_mask is None:
            var_mask = slice(None)
        elif self.data.isbacked and var_mask.dtype == bool:
            var_mask = var_mask.nonzero()[0]
        X = self.data.X[obs_mask, var_mask]
        return X

    def get_X_approximate_distribution(self) -> XApproximateDistribution:
        """return the approximate distribution of the X matrix."""
        if self.X_approximate_distribution is None:
            """Not yet evaluated."""
            assert self.dataset_config.X_approximate_distribution == "auto"
            self.data = self.data.to_memory()  # loads data
            self.X_approximate_distribution = estimate_distribution.estimate_approximate_distribution(self.data.X)

        return self.X_approximate_distribution

    def get_shape(self):
        return self.data.shape

    def query_var_array(self, term_name):
        return getattr(self.data.var, term_name)

    def query_obs_array(self, term_name):
        return getattr(self.data.obs, term_name)

    def get_obs_index(self):
        name = self.server_config.single_dataset__obs_names
        if name is None:
            return self.original_obs_index
        else:
            return self.data.obs[name]

    def get_obs_columns(self):
        return self.data.obs.columns

    def get_obs_keys(self):
        # return list of keys
        return self.data.obs.keys().to_list()

    def get_var_keys(self):
        # return list of keys
        return self.data.var.keys().to_list()

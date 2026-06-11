"""Utilities for reclustering selected observations.

This module intentionally builds small temporary AnnData objects instead of
copying the full selected expression matrix whenever possible.

Two compute paths are supported:

* representation mode: use an existing ``adata.obsm`` representation such as
  ``X_pca``. This is the most memory-friendly path.
* expression-gene mode: use all genes or resolve a whitelist/blacklist, slice
  the selected cells and genes, recompute PCA on that expression matrix, then
  run neighbors/Leiden/UMAP. This is required when the requested feature set
  should affect clustering.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fnmatch import fnmatchcase
import re
import time
from typing import Callable, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


OUTSIDE_RECLUSTER_CATEGORY = "__outside_recluster__"
GENE_FILTER_NONE = "none"
GENE_FILTER_ALL = "all"
GENE_FILTER_WHITELIST = "whitelist"
GENE_FILTER_BLACKLIST = "blacklist"
GENE_FILTER_MODES = {
    GENE_FILTER_NONE,
    GENE_FILTER_ALL,
    GENE_FILTER_WHITELIST,
    GENE_FILTER_BLACKLIST,
}


@dataclass(frozen=True)
class ReclusterResult:
    """A lightweight materialized reclustering result.

    The parent AnnData is not mutated. ``obs_indices`` maps rows in the
    reclustered result back to rows in the full, original AnnData. ``var_indices``
    is set only for expression-gene reclustering and maps the features used for PCA
    back to rows in the original ``var`` axis.
    """

    result_id: str
    user_id: str
    obs_indices: np.ndarray
    embedding: np.ndarray
    leiden: np.ndarray
    categories: list[str]
    params: dict
    var_indices: Optional[np.ndarray] = None
    var_names: Optional[list[str]] = None
    created_at: float = field(default_factory=time.time)

    @property
    def layout_name(self) -> str:
        return f"recluster_umap_{self.result_id}"

    @property
    def cluster_name(self) -> str:
        return f"recluster_leiden_{self.result_id}"

    def schema_payload(self) -> dict:
        """Return frontend schema for this generated result."""

        return {
            "result_id": self.result_id,
            "n_obs": int(self.obs_indices.size),
            "n_vars": None if self.var_indices is None else int(self.var_indices.size),
            "layout": {
                "name": self.layout_name,
                "type": "float32",
                "dims": [f"{self.layout_name}_0", f"{self.layout_name}_1"],
                "reclusterResultId": self.result_id,
                "reclusterField": "layout",
            },
            "cluster": {
                "name": self.cluster_name,
                "type": "categorical",
                "writable": False,
                "categories": self.categories + [OUTSIDE_RECLUSTER_CATEGORY],
                "reclusterResultId": self.result_id,
                "reclusterField": "leiden",
            },
            "params": self.params,
            "export_h5ad": f"recluster/obs/results/{self.result_id}/h5ad",
        }


@dataclass
class ReclusterJob:
    job_id: str
    user_id: str
    n_obs: int
    n_vars: Optional[int] = None
    gene_filter: Optional[dict] = None
    status: str = "queued"
    stage: str = "Queued"
    progress: float = 0.0
    result: Optional[dict] = None
    error: Optional[str] = None
    created_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)

    def update(self, *, status: Optional[str] = None, stage: Optional[str] = None, progress: Optional[float] = None):
        if status is not None:
            self.status = status
        if stage is not None:
            self.stage = stage
        if progress is not None:
            self.progress = float(progress)
        self.updated_at = time.time()

    def to_dict(self) -> dict:
        payload = {
            "job_id": self.job_id,
            "status": self.status,
            "stage": self.stage,
            "progress": self.progress,
            "n_obs": self.n_obs,
        }
        if self.n_vars is not None:
            payload["n_vars"] = self.n_vars
        if self.gene_filter is not None:
            payload["gene_filter"] = self.gene_filter
        if self.result is not None:
            payload["result"] = self.result
        if self.error is not None:
            payload["error"] = self.error
        return payload


ProgressCallback = Callable[[str, float], None]


def _report(progress: Optional[ProgressCallback], stage: str, value: float):
    if progress is not None:
        progress(stage, value)


def parse_gene_terms(value) -> list[str]:
    """Parse a UI-provided gene list into unique terms, preserving order.

    Accepts a Python list/tuple or a string containing newlines, commas,
    semicolons, tabs, or spaces. Empty terms are ignored. Shell-style wildcards
    such as ``MT-*`` and ``RPL*`` are supported later by ``resolve_gene_filter``.
    """

    if value is None:
        raw_terms = []
    elif isinstance(value, str):
        raw_terms = re.split(r"[\s,;]+", value)
    elif isinstance(value, (list, tuple, set)):
        raw_terms = list(value)
    else:
        raise ValueError("gene_list must be a string or a list of gene names")

    terms: list[str] = []
    seen = set()
    for term in raw_terms:
        term = str(term).strip()
        if not term or term in seen:
            continue
        terms.append(term)
        seen.add(term)
    return terms


def _normalise_match_value(value: str, *, case_sensitive: bool) -> str:
    return value if case_sensitive else value.upper()


def _is_glob_term(term: str) -> bool:
    return any(ch in term for ch in "*?[")


def resolve_gene_filter(var_names, params: Optional[dict]) -> Optional[dict]:
    """Resolve whitelist/blacklist terms against frontend-visible gene names.

    Returns ``None`` for no filter. Otherwise returns a serialisable summary plus
    a private ``var_indices`` key consumed by the backend. The private key is
    removed before returning summaries to the UI.
    """

    params = dict(params or {})
    mode = str(params.get("gene_filter_mode", GENE_FILTER_NONE) or GENE_FILTER_NONE).lower()
    if mode not in GENE_FILTER_MODES:
        raise ValueError("gene_filter_mode must be one of: none, all, whitelist, blacklist")

    terms = parse_gene_terms(params.get("gene_list", []))
    if mode == GENE_FILTER_NONE:
        if terms:
            raise ValueError("gene_list was provided but gene_filter_mode is 'none'")
        return None

    case_sensitive = bool(params.get("gene_filter_case_sensitive", False))
    names = np.asarray([str(name) for name in var_names], dtype=object)

    if mode == GENE_FILTER_ALL:
        if names.size < 2:
            raise ValueError("At least two genes are required for all-gene reclustering")
        return {
            "mode": mode,
            "terms": [],
            "case_sensitive": case_sensitive,
            "matched_gene_count": int(names.size),
            "included_gene_count": int(names.size),
            "unmatched_terms": [],
            "term_match_counts": {},
            "var_indices": np.arange(names.size, dtype=np.int64),
            "included_gene_names": [str(name) for name in names],
        }

    if not terms:
        raise ValueError("gene_list must contain at least one gene or wildcard term")

    match_names = np.asarray([_normalise_match_value(name, case_sensitive=case_sensitive) for name in names])

    matched_mask = np.zeros(names.shape[0], dtype=bool)
    term_match_counts: dict[str, int] = {}
    unmatched_terms: list[str] = []

    for term in terms:
        match_term = _normalise_match_value(term, case_sensitive=case_sensitive)
        if _is_glob_term(match_term):
            term_mask = np.fromiter((fnmatchcase(name, match_term) for name in match_names), dtype=bool, count=names.size)
        else:
            term_mask = match_names == match_term
        count = int(term_mask.sum())
        term_match_counts[term] = count
        if count == 0:
            unmatched_terms.append(term)
        matched_mask |= term_mask

    matched_indices = np.flatnonzero(matched_mask).astype(np.int64)
    if mode == GENE_FILTER_WHITELIST:
        var_indices = matched_indices
    else:
        var_indices = np.flatnonzero(~matched_mask).astype(np.int64)

    if var_indices.size < 2:
        raise ValueError(
            "Gene filter leaves fewer than two genes available for reclustering. "
            "Use a larger whitelist or a smaller blacklist."
        )

    return {
        "mode": mode,
        "terms": terms,
        "case_sensitive": case_sensitive,
        "matched_gene_count": int(matched_indices.size),
        "included_gene_count": int(var_indices.size),
        "unmatched_terms": unmatched_terms,
        "term_match_counts": term_match_counts,
        "var_indices": var_indices,
        "included_gene_names": [str(names[i]) for i in var_indices],
    }


def public_gene_filter_summary(gene_filter: Optional[dict]) -> Optional[dict]:
    """Return a JSON-safe gene-selection summary without private arrays."""

    if gene_filter is None:
        return None
    return {
        key: value
        for key, value in gene_filter.items()
        if key not in {"var_indices", "included_gene_names"}
    }


def normalize_recluster_params(params: Optional[dict], *, n_obs: int, rep_dim: int) -> dict:
    """Validate and normalize representation-based reclustering parameters."""

    params = dict(params or {})
    use_rep = str(params.get("use_rep", "X_pca"))

    n_neighbors = int(params.get("n_neighbors", 15))
    if n_neighbors < 2:
        raise ValueError("n_neighbors must be at least 2")
    if n_neighbors >= n_obs:
        n_neighbors = n_obs - 1

    n_pcs = params.get("n_pcs", 50)
    if n_pcs is None or n_pcs == "":
        n_pcs = rep_dim
    n_pcs = int(n_pcs)
    if n_pcs < 1:
        raise ValueError("n_pcs must be at least 1")
    n_pcs = min(n_pcs, rep_dim)

    resolution = float(params.get("resolution", 1.0))
    if resolution <= 0:
        raise ValueError("resolution must be greater than 0")

    min_dist = float(params.get("min_dist", 0.5))
    if min_dist < 0:
        raise ValueError("min_dist must be non-negative")

    random_state = int(params.get("random_state", 0))

    return {
        "compute_mode": "obsm",
        "use_rep": use_rep,
        "n_neighbors": n_neighbors,
        "n_pcs": n_pcs,
        "resolution": resolution,
        "min_dist": min_dist,
        "random_state": random_state,
        "gene_filter_mode": GENE_FILTER_NONE,
    }


def normalize_recluster_gene_params(params: Optional[dict], *, n_obs: int, n_vars: int) -> dict:
    """Validate and normalize expression-gene reclustering parameters."""

    params = dict(params or {})
    n_neighbors = int(params.get("n_neighbors", 15))
    if n_neighbors < 2:
        raise ValueError("n_neighbors must be at least 2")
    if n_neighbors >= n_obs:
        n_neighbors = n_obs - 1

    n_pcs = params.get("n_pcs", 50)
    if n_pcs is None or n_pcs == "":
        n_pcs = min(n_obs - 1, n_vars)
    n_pcs = int(n_pcs)
    if n_pcs < 1:
        raise ValueError("n_pcs must be at least 1")

    resolution = float(params.get("resolution", 1.0))
    if resolution <= 0:
        raise ValueError("resolution must be greater than 0")

    min_dist = float(params.get("min_dist", 0.5))
    if min_dist < 0:
        raise ValueError("min_dist must be non-negative")

    random_state = int(params.get("random_state", 0))
    log1p = bool(params.get("gene_filter_log1p", False))
    scale = bool(params.get("gene_filter_scale", False))

    # svd_solver='arpack' requires n_comps < min(n_obs, n_vars). Cap one below
    # that bound. If the selected feature set is tiny, we skip PCA and use X.
    max_pca_comps = max(1, min(n_obs, n_vars) - 1)
    n_pcs = min(n_pcs, max_pca_comps)

    mode = str(params.get("gene_filter_mode", GENE_FILTER_WHITELIST) or GENE_FILTER_WHITELIST).lower()
    return {
        "compute_mode": "expression_genes",
        "gene_filter_mode": mode,
        "n_neighbors": n_neighbors,
        "n_pcs": n_pcs,
        "resolution": resolution,
        "min_dist": min_dist,
        "random_state": random_state,
        "gene_filter_log1p": log1p,
        "gene_filter_scale": scale,
    }


def _slice_expression_matrix(X, obs_indices: np.ndarray, var_indices: np.ndarray):
    """Slice X to selected observations and variables.

    This function avoids ``adata.copy()``. Sparse inputs remain sparse. Dense and
    backed dense inputs are materialized as a selected-cell by selected-gene
    float32 array, which is unavoidable if the clustering feature set is changed.
    """

    if sparse.issparse(X):
        return X[obs_indices, :][:, var_indices].astype(np.float32)

    # For regular ndarray this gives the smallest direct dense slice. Some backed
    # array implementations do not support two-axis fancy indexing; fall back to
    # chunking by selected rows in that case.
    try:
        sliced = X[np.ix_(obs_indices, var_indices)]
        if sparse.issparse(sliced):
            return sliced.astype(np.float32)
        return np.asarray(sliced, dtype=np.float32, order="C")
    except Exception:
        chunks = []
        chunk_size = 1024
        for start in range(0, obs_indices.size, chunk_size):
            obs_chunk = obs_indices[start : start + chunk_size]
            try:
                chunk = X[np.ix_(obs_chunk, var_indices)]
            except Exception:
                chunk = X[obs_chunk, :]
                chunk = chunk[:, var_indices]
            chunks.append(chunk)

        if any(sparse.issparse(chunk) for chunk in chunks):
            return sparse.vstack([chunk if sparse.issparse(chunk) else sparse.csr_matrix(chunk) for chunk in chunks]).astype(
                np.float32
            )
        return np.asarray(np.vstack(chunks), dtype=np.float32, order="C")


def _build_obs_frame(source_adata, obs_indices: np.ndarray) -> pd.DataFrame:
    obs_names = source_adata.obs_names.to_numpy()[obs_indices]
    obs = pd.DataFrame(index=pd.Index([str(x) for x in obs_names], name="obs_names"))
    obs["cellxgene_original_obs_index"] = obs_indices
    return obs


def _run_graph_clustering(
    work: ad.AnnData,
    *,
    result_id: str,
    params: dict,
    use_rep: str,
    neighbors_n_pcs: Optional[int],
    progress: Optional[ProgressCallback],
) -> tuple[np.ndarray, np.ndarray, list[str], dict]:
    """Run neighbors -> Leiden -> UMAP on an already prepared work AnnData."""

    neighbors_key = f"cxg_recluster_neighbors_{result_id}"
    leiden_key = f"cxg_recluster_leiden_{result_id}"
    umap_key = f"X_cxg_recluster_umap_{result_id}"

    _report(progress, "Computing neighbors", 0.25)
    sc.pp.neighbors(
        work,
        n_neighbors=params["n_neighbors"],
        n_pcs=neighbors_n_pcs,
        use_rep=use_rep,
        key_added=neighbors_key,
        random_state=params["random_state"],
        copy=False,
    )

    _report(progress, "Computing Leiden clusters", 0.55)
    sc.tl.leiden(
        work,
        resolution=params["resolution"],
        key_added=leiden_key,
        neighbors_key=neighbors_key,
        random_state=params["random_state"],
        flavor="igraph",
        n_iterations=2,
        directed=False,
        copy=False,
    )

    _report(progress, "Computing UMAP", 0.75)
    sc.tl.umap(
        work,
        neighbors_key=neighbors_key,
        key_added=umap_key,
        min_dist=params["min_dist"],
        random_state=params["random_state"],
        copy=False,
    )

    _report(progress, "Registering result", 0.95)
    leiden_series = work.obs[leiden_key].astype("category")
    categories = [str(x) for x in leiden_series.cat.categories.to_list()]

    extra_params = {
        "neighbors_key": neighbors_key,
        "leiden_key": leiden_key,
        "umap_key": umap_key,
        "leiden_flavor": "igraph",
        "leiden_n_iterations": 2,
        "leiden_directed": False,
    }
    return (
        work.obsm[umap_key].astype(np.float32, copy=False),
        leiden_series.astype(str).to_numpy(),
        categories,
        extra_params,
    )


def recluster_from_obsm(
    source_adata,
    selected_obs_indices,
    *,
    result_id: str,
    user_id: str,
    params: Optional[dict] = None,
    progress: Optional[ProgressCallback] = None,
) -> ReclusterResult:
    """Recluster selected cells using an existing representation in ``.obsm``.

    The compute path stores only a selected-cell representation in memory, not a
    full selected copy of ``X``. The returned result is lightweight and can be
    served back to the frontend as generated layout/annotation columns.
    """

    obs_indices = np.asarray(selected_obs_indices, dtype=np.int64)
    obs_indices = np.unique(obs_indices)
    n_obs = int(obs_indices.size)
    if n_obs < 10:
        raise ValueError("Select at least 10 cells to recluster")

    params = dict(params or {})
    use_rep = str(params.get("use_rep", "X_pca"))
    if use_rep not in source_adata.obsm:
        available = ", ".join(sorted(str(k) for k in source_adata.obsm.keys()))
        raise ValueError(f"Representation {use_rep!r} was not found in adata.obsm. Available: {available}")

    rep = source_adata.obsm[use_rep]
    if len(rep.shape) != 2:
        raise ValueError(f"Representation {use_rep!r} must be 2-dimensional")
    rep_dim = int(rep.shape[1])
    params = normalize_recluster_params(params, n_obs=n_obs, rep_dim=rep_dim)

    _report(progress, "Preparing selected cells", 0.10)
    # Use only the requested latent dimensions. np.asarray(..., order='C')
    # materializes a small dense matrix needed by pynndescent/umap.
    X_rep = rep[obs_indices, :]
    X_rep = np.asarray(X_rep[:, : params["n_pcs"]], dtype=np.float32, order="C")

    obs = _build_obs_frame(source_adata, obs_indices)

    # This is a working object whose X is the PCA/latent representation, not
    # the full expression matrix. It is intentionally discarded after extracting
    # the small result arrays.
    work = ad.AnnData(X=X_rep, obs=obs)

    embedding, leiden, categories, key_params = _run_graph_clustering(
        work,
        result_id=result_id,
        params=params,
        use_rep="X",
        neighbors_n_pcs=None,
        progress=progress,
    )

    return ReclusterResult(
        result_id=result_id,
        user_id=user_id,
        obs_indices=obs_indices,
        embedding=embedding,
        leiden=leiden,
        categories=categories,
        params={**params, **key_params},
    )


def recluster_from_expression_genes(
    source_adata,
    selected_obs_indices,
    *,
    var_indices,
    var_names,
    gene_filter: dict,
    result_id: str,
    user_id: str,
    params: Optional[dict] = None,
    progress: Optional[ProgressCallback] = None,
) -> ReclusterResult:
    """Recluster selected cells using expression genes.

    This path is used for all-gene, whitelist, and blacklist reclustering. It
    recomputes PCA on the selected expression matrix because existing ``X_pca``
    cannot be adjusted after the fact to remove, include, or reweight specific
    genes.
    """

    obs_indices = np.asarray(selected_obs_indices, dtype=np.int64)
    obs_indices = np.unique(obs_indices)
    var_indices = np.asarray(var_indices, dtype=np.int64)
    n_obs = int(obs_indices.size)
    n_vars = int(var_indices.size)
    if n_obs < 10:
        raise ValueError("Select at least 10 cells to recluster")
    if n_vars < 2:
        raise ValueError("At least two genes are required for expression-gene reclustering")

    params = normalize_recluster_gene_params(params, n_obs=n_obs, n_vars=n_vars)

    _report(progress, "Preparing selected genes", 0.10)
    X = _slice_expression_matrix(source_adata.X, obs_indices, var_indices)
    obs = _build_obs_frame(source_adata, obs_indices)
    var = pd.DataFrame(index=pd.Index([str(x) for x in var_names], name="var_names"))
    var["cellxgene_original_var_index"] = var_indices
    work = ad.AnnData(X=X, obs=obs, var=var)

    if params["gene_filter_log1p"]:
        _report(progress, "Applying log1p transform", 0.16)
        sc.pp.log1p(work, copy=False)

    if params["gene_filter_scale"]:
        _report(progress, "Scaling selected genes", 0.19)
        # zero_center=False keeps sparse matrices sparse and avoids a large dense
        # copy. It is therefore more appropriate for an interactive workflow.
        sc.pp.scale(work, zero_center=False, copy=False)

    pca_key = None
    use_rep = "X"
    neighbors_n_pcs = None
    if params["n_pcs"] >= 2 and n_vars >= 3 and n_obs >= 3:
        _report(progress, "Computing PCA from selected genes", 0.22)
        sc.pp.pca(
            work,
            n_comps=params["n_pcs"],
            random_state=params["random_state"],
            svd_solver="arpack",
            copy=False,
        )
        pca_key = "X_pca"
        use_rep = "X_pca"
        neighbors_n_pcs = params["n_pcs"]

    embedding, leiden, categories, key_params = _run_graph_clustering(
        work,
        result_id=result_id,
        params=params,
        use_rep=use_rep,
        neighbors_n_pcs=neighbors_n_pcs,
        progress=progress,
    )

    gene_filter_summary = public_gene_filter_summary(gene_filter)
    result_params = {
        **params,
        **key_params,
        "pca_key": pca_key,
        "gene_filter": gene_filter_summary,
    }

    return ReclusterResult(
        result_id=result_id,
        user_id=user_id,
        obs_indices=obs_indices,
        var_indices=var_indices,
        var_names=[str(x) for x in var_names],
        embedding=embedding,
        leiden=leiden,
        categories=categories,
        params=result_params,
    )

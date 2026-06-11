import base64
import csv
import fnmatch
import io
import json
from contextlib import contextmanager

import numpy as np
import pandas as pd
from scipy import sparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402


PLOT_SELECTION_OBS = "cellxgene_plot_selection"
DE_GROUP_OBS = "cellxgene_de_group"

_IMAGE_MIMETYPES = {
    "png": "image/png",
    "svg": "image/svg+xml",
    "pdf": "application/pdf",
}

_BLOCKED_KWARGS = {
    "adata",
    "show",
    "save",
    "ax",
    "return_fig",
}


SCANPY_PLOT_DEFINITIONS = [
    {
        "id": "embedding",
        "label": "Embedding / UMAP / t-SNE / PCA",
        "category": "embedding",
        "description": "Scatter plot on an embedding using scanpy.pl.embedding.",
    },
    {"id": "dotplot", "label": "Dot plot", "category": "gene_by_group"},
    {"id": "matrixplot", "label": "Matrix plot", "category": "gene_by_group"},
    {"id": "heatmap", "label": "Heatmap", "category": "gene_by_group"},
    {"id": "violin", "label": "Violin", "category": "gene_by_group"},
    {"id": "stacked_violin", "label": "Stacked violin", "category": "gene_by_group"},
    {"id": "tracksplot", "label": "Tracks plot", "category": "gene_by_group"},
    {"id": "embedding_density", "label": "Embedding density", "category": "embedding"},
    {"id": "density_scatter", "label": "Density scatter", "category": "custom"},
    {"id": "stacked_barplot", "label": "Stacked barplot", "category": "custom"},
    {"id": "highest_expr_genes", "label": "Highest expressed genes", "category": "qc"},
    {"id": "highly_variable_genes", "label": "Highly variable genes", "category": "qc"},
    {"id": "pca_variance_ratio", "label": "PCA variance ratio", "category": "tools"},
    {"id": "paga", "label": "PAGA", "category": "tools"},
    {"id": "dendrogram", "label": "Dendrogram", "category": "tools"},
]

DE_PLOT_DEFINITIONS = [
    {"id": "rank_genes_groups", "label": "Rank genes groups"},
    {"id": "rank_genes_groups_dotplot", "label": "Rank genes groups dot plot"},
    {"id": "rank_genes_groups_matrixplot", "label": "Rank genes groups matrix plot"},
    {"id": "rank_genes_groups_stacked_violin", "label": "Rank genes groups stacked violin"},
    {"id": "rank_genes_groups_tracksplot", "label": "Rank genes groups tracks plot"},
    {"id": "rank_genes_groups_heatmap", "label": "Rank genes groups heatmap"},
    {"id": "rank_genes_groups_violin", "label": "Rank genes groups violin"},
    {"id": "volcano", "label": "Volcano plot"},
    {"id": "de_volcano", "label": "Volcano plot"},
]

DE_GROUPED_PLOT_IDS = {
    "rank_genes_groups_dotplot",
    "rank_genes_groups_matrixplot",
    "rank_genes_groups_stacked_violin",
    "rank_genes_groups_tracksplot",
    "rank_genes_groups_heatmap",
}


class ScanpyPlotError(ValueError):
    pass


def parse_gene_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        terms = value
    else:
        for sep in [",", ";", "\t", "\n", "\r"]:
            value = str(value).replace(sep, " ")
        terms = value.split(" ")
    return [str(term).strip() for term in terms if str(term).strip()]


def _term_matches_gene_name(gene_name, term, *, case_sensitive=False):
    gene_name = str(gene_name)
    term = str(term)
    if not case_sensitive:
        gene_name = gene_name.lower()
        term = term.lower()
    if any(char in term for char in ("*", "?", "[")):
        return fnmatch.fnmatchcase(gene_name, term)
    return gene_name == term


def resolve_de_gene_exclusion(var_names, value, *, case_sensitive=False):
    terms = parse_gene_list(value)
    names = np.asarray([str(name) for name in var_names], dtype=object)
    excluded_mask = np.zeros(names.shape[0], dtype=bool)
    unmatched_terms = []

    for term in terms:
        term_mask = np.asarray(
            [
                _term_matches_gene_name(
                    gene_name,
                    term,
                    case_sensitive=case_sensitive,
                )
                for gene_name in names
            ],
            dtype=bool,
        )
        if not np.any(term_mask):
            unmatched_terms.append(term)
        excluded_mask = np.logical_or(excluded_mask, term_mask)

    return {
        "terms": terms,
        "excluded_mask": excluded_mask,
        "excluded_count": int(np.count_nonzero(excluded_mask)),
        "unmatched_terms": unmatched_terms,
    }


def apply_de_gene_exclusion(adata, settings):
    exclusion = resolve_de_gene_exclusion(
        adata.var_names,
        settings.get("de_exclude_genes"),
        case_sensitive=_as_bool(settings.get("de_exclude_case_sensitive"), False),
    )
    if not exclusion["terms"]:
        return {
            "terms": [],
            "excluded_count": 0,
            "remaining_gene_count": int(adata.n_vars),
            "unmatched_terms": [],
        }

    if exclusion["excluded_count"] == 0:
        preview = ", ".join(exclusion["terms"][:10])
        raise ScanpyPlotError(
            "DE exclude gene list did not match any genes. "
            f"Check gene symbols or wildcard patterns: {preview}"
        )

    include_mask = np.logical_not(exclusion["excluded_mask"])
    remaining_gene_count = int(np.count_nonzero(include_mask))
    if remaining_gene_count < 2:
        raise ScanpyPlotError(
            "DE gene exclusion leaves fewer than two genes available. "
            "Use fewer excluded genes or remove overly broad wildcards."
        )

    # If use_raw is automatic, Scanpy may rank against adata.raw and bypass the
    # var subset below. Force the temporary X/layer matrix so the blacklist is
    # respected for this DE run.
    settings["_de_force_use_raw_false"] = True
    adata._inplace_subset_var(include_mask)

    return {
        "terms": exclusion["terms"],
        "excluded_count": int(exclusion["excluded_count"]),
        "remaining_gene_count": remaining_gene_count,
        "unmatched_terms": exclusion["unmatched_terms"],
    }


def parse_group_list(value):
    """Parse a group list without splitting on spaces inside category names."""
    if value is None:
        return []
    if isinstance(value, list):
        terms = value
    else:
        value = str(value)
        for sep in [";", "\t", "\n", "\r"]:
            value = value.replace(sep, ",")
        terms = value.split(",")
    return [str(term).strip() for term in terms if str(term).strip()]


def _json_kwargs(value):
    if not value:
        return {}
    if isinstance(value, dict):
        kwargs = value
    else:
        kwargs = json.loads(value)
    if not isinstance(kwargs, dict):
        raise ScanpyPlotError("Advanced Scanpy kwargs must be a JSON object.")
    return {str(k): v for k, v in kwargs.items() if str(k) not in _BLOCKED_KWARGS}


def _as_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes", "on")
    return bool(value)


def _as_optional_float(value):
    if value is None or value == "":
        return None
    return float(value)


def _as_optional_int(value):
    if value is None or value == "":
        return None
    return int(value)


def _as_optional_bool_select(value):
    if value is None or value == "" or value == "auto":
        return None
    return _as_bool(value)


def _as_optional_categories(value):
    categories = parse_group_list(value)
    return categories or None


def _figsize(settings):
    return (float(settings.get("width") or 7), float(settings.get("height") or 5))


def _title(settings):
    return settings.get("title") or None


def _listify_color(value):
    if value is None:
        return None
    if isinstance(value, list):
        values = [str(item).strip() for item in value if str(item).strip()]
        return values or None
    values = parse_gene_list(value)
    return values or None



def _get_vector(adata, name, *, layer=None, use_raw=None):
    """Return a dense 1-D vector for a gene/var name or obs annotation."""
    name = str(name)
    if name in adata.obs.columns:
        values = adata.obs[name]
        if isinstance(values.dtype, pd.CategoricalDtype):
            return values.astype(str).to_numpy()
        return np.asarray(values)

    if name not in set(map(str, adata.var_names)):
        raise ScanpyPlotError(f"Unknown gene or annotation: {name}")

    source = adata.raw if use_raw and adata.raw is not None else adata
    var_names = np.asarray([str(v) for v in source.var_names], dtype=object)
    matches = np.flatnonzero(var_names == name)
    if matches.size == 0:
        raise ScanpyPlotError(f"Gene {name!r} is not available in the selected data source.")
    idx = int(matches[0])

    if layer and use_raw:
        raise ScanpyPlotError("Layer and use_raw cannot both be used for custom plots.")
    matrix = source.X if not layer else adata.layers[layer]
    values = matrix[:, idx]
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).reshape(-1)


def _as_numeric_vector(values, label):
    try:
        numeric = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    except Exception as exc:
        raise ScanpyPlotError(f"{label} must be numeric for this plot.") from exc
    if not np.isfinite(numeric).any():
        raise ScanpyPlotError(f"{label} has no finite numeric values for this plot.")
    return numeric


def _categorical_series_from_obs_or_gene(adata, name, *, layer=None, use_raw=None, cutoff=0.0):
    name = str(name or "").strip()
    if not name:
        raise ScanpyPlotError("Choose an annotation or gene.")
    if name in adata.obs.columns:
        series = adata.obs[name]
        if not isinstance(series.dtype, pd.CategoricalDtype):
            series = series.astype("category")
        return series.cat.remove_unused_categories()

    values = _as_numeric_vector(
        _get_vector(adata, name, layer=layer, use_raw=use_raw),
        name,
    )
    labels = np.where(values > float(cutoff), f"{name} > {cutoff:g}", f"{name} <= {cutoff:g}")
    return pd.Series(pd.Categorical(labels), index=adata.obs_names)


def _ordered_categories(series, order_text=None):
    categories = [str(c) for c in series.cat.categories]
    requested = parse_group_list(order_text)
    if requested:
        missing = [value for value in requested if value not in categories]
        if missing:
            raise ScanpyPlotError(f"Unknown category in order list: {', '.join(missing)}")
        categories = requested + [value for value in categories if value not in requested]
    return categories


def _first_rank_group(adata, de_info=None, requested_group=None):
    requested_group = str(requested_group or "").strip()
    names = adata.uns.get("rank_genes_groups", {}).get("names")
    available = list(names.dtype.names or []) if hasattr(names, "dtype") else []
    if requested_group:
        if requested_group not in available:
            raise ScanpyPlotError(f"Unknown DE volcano group: {requested_group}")
        return requested_group
    csv_group = (de_info or {}).get("csv_group")
    if csv_group:
        return csv_group
    groups_to_plot = (de_info or {}).get("groups_to_plot")
    if groups_to_plot:
        return groups_to_plot[0]
    if available:
        return available[0]
    return None


def _rank_genes_groups_df(adata, group=None):
    try:
        return sc.get.rank_genes_groups_df(adata, group=group)
    except Exception:
        return sc.get.rank_genes_groups_df(adata, group=None)


def _json_safe_value(value):
    if value is None or pd.isna(value):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _dataframe_csv_and_preview(df, max_preview_rows=25):
    export_df = df.copy()
    buf = io.StringIO()
    export_df.to_csv(buf, index=False, quoting=csv.QUOTE_MINIMAL)
    preview_df = export_df.head(max_preview_rows).replace([np.inf, -np.inf], np.nan)
    columns = [str(column) for column in preview_df.columns]
    rows = []
    for record in preview_df.to_dict(orient="records"):
        rows.append({str(key): _json_safe_value(value) for key, value in record.items()})
    return buf.getvalue(), {
        "columns": columns,
        "rows": rows,
        "row_count": int(export_df.shape[0]),
    }

def _validate_genes(adata, genes, required=False):
    genes = parse_gene_list(genes)
    if required and not genes:
        raise ScanpyPlotError("Enter at least one gene for this plot.")
    if not genes:
        return []

    var_names = set(map(str, adata.var_names))
    obs_names = set(map(str, adata.obs.columns))
    missing = [gene for gene in genes if gene not in var_names and gene not in obs_names]
    if missing:
        preview = ", ".join(missing[:10])
        extra = "" if len(missing) <= 10 else f", and {len(missing) - 10} more"
        raise ScanpyPlotError(f"Unknown gene or annotation for plotting: {preview}{extra}")
    return genes


def _validate_groupby(adata, groupby, required=False):
    if not groupby:
        if required:
            raise ScanpyPlotError("Choose a categorical annotation for groupby.")
        return None
    if groupby not in adata.obs.columns:
        raise ScanpyPlotError(f"Unknown obs annotation for groupby: {groupby}")
    return groupby


def _ensure_categorical_without_unused(adata, column):
    if not column or column not in adata.obs:
        return []
    if not isinstance(adata.obs[column].dtype, pd.CategoricalDtype):
        adata.obs[column] = adata.obs[column].astype("category")
    adata.obs[column] = adata.obs[column].cat.remove_unused_categories()
    return [str(category) for category in adata.obs[column].cat.categories]


def _validate_basis(adata, basis):
    if not basis:
        raise ScanpyPlotError("Choose an embedding basis.")
    key = f"X_{basis}"
    if key not in adata.obsm:
        available = sorted(k[2:] for k in adata.obsm.keys() if isinstance(k, str) and k.startswith("X_"))
        raise ScanpyPlotError(f"Embedding {basis!r} is not available. Available embeddings: {', '.join(available)}")
    return basis


def _valid_export_format(fmt):
    fmt = (fmt or "png").lower().strip().lstrip(".")
    if fmt not in _IMAGE_MIMETYPES:
        raise ScanpyPlotError("Plot export format must be png, svg, or pdf.")
    return fmt


def _prepare_x_for_plotting(adata, settings):
    layer = settings.get("layer") or None
    if layer and layer not in adata.layers:
        raise ScanpyPlotError(f"Layer {layer!r} is not available in this dataset.")

    if _as_bool(settings.get("log1p"), False):
        sc.pp.log1p(adata, layer=layer)

    return layer


def _common_plot_scale_kwargs(settings):
    kwargs = {
        "vmin": _as_optional_float(settings.get("vmin")),
        "vmax": _as_optional_float(settings.get("vmax")),
        "vcenter": _as_optional_float(settings.get("vcenter")),
    }
    return {k: v for k, v in kwargs.items() if v is not None}


def _common_grouped_plot_kwargs(settings):
    kwargs = {
        "standard_scale": settings.get("standard_scale") or None,
        "swap_axes": _as_bool(settings.get("swap_axes"), False),
        "num_categories": _as_optional_int(settings.get("num_categories")),
        "categories_order": _as_optional_categories(settings.get("categories_order")),
        "figsize": _figsize(settings),
        "title": _title(settings),
        "show": False,
    }
    kwargs.update(_common_plot_scale_kwargs(settings))
    return {k: v for k, v in kwargs.items() if v is not None}


@contextmanager
def _figure_context(settings):
    plt.close("all")
    dpi = int(settings.get("dpi") or 150)
    width = float(settings.get("width") or 7)
    height = float(settings.get("height") or 5)
    font_size = float(settings.get("font_size") or 10)
    with plt.rc_context({"figure.dpi": dpi, "savefig.dpi": dpi, "font.size": font_size}):
        fig = plt.figure(figsize=(width, height), dpi=dpi)
        try:
            yield fig
        finally:
            pass


def _figure_from_return(ret):
    if ret is None:
        return plt.gcf()
    if hasattr(ret, "make_figure"):
        try:
            ret.make_figure()
        except Exception:
            pass
        if hasattr(ret, "fig") and ret.fig is not None:
            return ret.fig
    if hasattr(ret, "fig") and ret.fig is not None:
        return ret.fig
    if hasattr(ret, "figure") and ret.figure is not None:
        return ret.figure
    if isinstance(ret, dict):
        for value in ret.values():
            fig = _figure_from_return(value)
            if fig is not None:
                return fig
    if isinstance(ret, (list, tuple)):
        for value in ret:
            fig = _figure_from_return(value)
            if fig is not None:
                return fig
    return plt.gcf()


def _save_figure(fig, fmt, settings):
    buf = io.BytesIO()
    kwargs = {
        "format": fmt,
        "dpi": int(settings.get("dpi") or 150),
        "bbox_inches": "tight",
        "transparent": _as_bool(settings.get("transparent"), False),
    }
    if fmt == "png":
        kwargs["facecolor"] = "white"
    fig.savefig(buf, **kwargs)
    return buf.getvalue()


def _build_response(fig, settings, filename_prefix, csv_text=None, csv_preview=None):
    export_format = _valid_export_format(settings.get("format"))
    preview_bytes = _save_figure(fig, "png", settings)
    download_bytes = preview_bytes if export_format == "png" else _save_figure(fig, export_format, settings)
    plt.close(fig)

    response = {
        "preview_mime": "image/png",
        "preview_base64": base64.b64encode(preview_bytes).decode("ascii"),
        "download_mime": _IMAGE_MIMETYPES[export_format],
        "download_base64": base64.b64encode(download_bytes).decode("ascii"),
        "download_filename": f"{filename_prefix}.{export_format}",
        "format": export_format,
    }
    if csv_text is not None:
        response["csv_mime"] = "text/csv"
        response["csv_base64"] = base64.b64encode(csv_text.encode("utf-8")).decode("ascii")
        response["csv_filename"] = f"{filename_prefix}.csv"
    if csv_preview is not None:
        response["csv_preview"] = csv_preview
    return response


def _csv_from_rank_genes_groups(adata, group=None):
    df = _rank_genes_groups_df(adata, group=group)
    return _dataframe_csv_and_preview(df)


def _set_figure_size(fig, settings):
    width = float(settings.get("width") or 7)
    height = float(settings.get("height") or 5)
    fig.set_size_inches(width, height, forward=True)



def _run_embedding_density_plot(adata, settings, common):
    basis = _validate_basis(adata, settings.get("basis"))
    groupby = _validate_groupby(
        adata,
        settings.get("density_groupby") or settings.get("embedding_density_groupby"),
        required=False,
    )
    if groupby:
        _ensure_categorical_without_unused(adata, groupby)
    key = (
        settings.get("density_key")
        or settings.get("embedding_density_key")
        or f"{basis}_density" + (f"_{groupby}" if groupby else "")
    )
    sc.tl.embedding_density(
        adata,
        basis=basis,
        groupby=groupby,
        key_added=key,
        components=settings.get("components") or None,
    )
    if groupby:
        group = settings.get("density_group") or settings.get("embedding_density_group") or "all"
        if str(group).lower() == "none":
            group = None
        elif isinstance(group, str) and "," in group:
            group = parse_group_list(group)
    else:
        # Overall density has no categorical covariate.  Scanpy's plotting API
        # accepts either key or groupby, not both; when using the explicit key
        # generated above, do not also pass groupby.  Without groupby, group
        # should be None to request the overall density instead of category panels.
        group = None
    kwargs = {
        "basis": basis,
        "key": key,
        "group": group,
        "color_map": settings.get("color_map") or "YlOrRd",
        "bg_dotsize": _as_optional_float(
            settings.get("density_bg_dotsize") or settings.get("embedding_density_bg_dotsize")
        ),
        "fg_dotsize": _as_optional_float(
            settings.get("density_fg_dotsize") or settings.get("embedding_density_fg_dotsize")
        ),
        "vmin": _as_optional_float(settings.get("vmin")),
        "vmax": _as_optional_float(settings.get("vmax")),
        "vcenter": _as_optional_float(settings.get("vcenter")),
        "ncols": _as_optional_int(settings.get("ncols")),
        "title": _title(settings),
        "show": False,
    }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    # The workbench controls these two values.  Allowing them through advanced
    # kwargs can recreate Scanpy's "either pass key or groupby but not both" error.
    common = {k: v for k, v in common.items() if k not in ("key", "groupby")}
    kwargs.update(common)
    return sc.pl.embedding_density(adata, **kwargs)


def _run_density_scatter_plot(adata, settings, *, layer=None, use_raw=None):
    genes = parse_gene_list(settings.get("density_scatter_genes") or settings.get("genes"))
    if len(genes) != 2:
        x_gene = str(settings.get("density_gene_x") or "").strip()
        y_gene = str(settings.get("density_gene_y") or "").strip()
        genes = [gene for gene in (x_gene, y_gene) if gene]
    if len(genes) != 2:
        raise ScanpyPlotError("Density scatter requires exactly two genes or numeric annotations.")
    x = _as_numeric_vector(_get_vector(adata, genes[0], layer=layer, use_raw=use_raw), genes[0])
    y = _as_numeric_vector(_get_vector(adata, genes[1], layer=layer, use_raw=use_raw), genes[1])
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        raise ScanpyPlotError("The two selected genes have no finite values to plot.")
    x = x[mask]
    y = y[mask]
    gridsize = int(settings.get("density_scatter_bins") or settings.get("density_hexagons") or 50)
    mincnt = int(settings.get("density_min_count") or 1)
    cutoff = _as_optional_float(
        settings.get("density_scatter_cutoff") or settings.get("density_expression_cutoff")
    )

    fig, ax = plt.subplots(figsize=_figsize(settings))
    hb = ax.hexbin(
        x,
        y,
        gridsize=gridsize,
        mincnt=mincnt,
        cmap=settings.get("color_map") or "viridis",
    )
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("count")
    if cutoff is not None:
        ax.axvline(cutoff, linestyle="--", linewidth=1)
        ax.axhline(cutoff, linestyle="--", linewidth=1)
    ax.set_xlabel(genes[0])
    ax.set_ylabel(genes[1])
    ax.set_title(_title(settings) or f"Density on {x.size} cells")
    ax.grid(True, alpha=0.25)
    return fig


def _run_stacked_barplot(adata, settings, *, layer=None, use_raw=None):
    x_name = settings.get("stacked_bar_x") or settings.get("groupby")
    color_name = settings.get("stacked_bar_color") or settings.get("stacked_bar_colorby") or settings.get("color")
    if not x_name or not color_name:
        raise ScanpyPlotError("Stacked barplot requires an x annotation/gene and a color annotation/gene.")
    cutoff = float(settings.get("stacked_bar_gene_cutoff") or settings.get("density_scatter_cutoff") or 0)
    x_series = _categorical_series_from_obs_or_gene(
        adata,
        x_name,
        layer=layer,
        use_raw=use_raw,
        cutoff=cutoff,
    )
    color_series = _categorical_series_from_obs_or_gene(
        adata,
        color_name,
        layer=layer,
        use_raw=use_raw,
        cutoff=cutoff,
    )
    x_order = _ordered_categories(x_series, settings.get("stacked_bar_x_order") or settings.get("categories_order"))
    color_order = _ordered_categories(color_series, settings.get("stacked_bar_color_order"))
    counts = pd.crosstab(x_series, color_series).reindex(index=x_order, columns=color_order, fill_value=0)
    mode = (settings.get("stacked_bar_mode") or "count").lower()

    fig, ax = plt.subplots(figsize=_figsize(settings))
    x_positions = np.arange(len(counts.index))
    if mode == "proportion":
        values = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        bottom = np.zeros(len(values.index))
        for category in values.columns:
            ax.bar(x_positions, values[category].to_numpy(), bottom=bottom, label=str(category))
            bottom += values[category].to_numpy()
        ax.set_ylabel("proportion")
    elif mode == "streamgraph":
        values = counts.T.to_numpy(dtype=float)
        ax.stackplot(x_positions, values, labels=[str(c) for c in counts.columns], baseline="sym")
        ax.set_ylabel("streamgraph count")
    else:
        values = counts
        bottom = np.zeros(len(values.index))
        for category in values.columns:
            ax.bar(x_positions, values[category].to_numpy(), bottom=bottom, label=str(category))
            bottom += values[category].to_numpy()
        ax.set_ylabel("count")

    rotation = float(settings.get("stacked_bar_label_rotation") or settings.get("rotation") or 60)
    font_size = _as_optional_float(settings.get("stacked_bar_x_font_size"))
    label_shift = _as_optional_float(settings.get("stacked_bar_label_shift"))
    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        [str(x) for x in counts.index],
        rotation=rotation,
        ha="right" if abs(rotation) > 0 else "center",
        fontsize=font_size,
    )
    if label_shift is not None:
        ax.tick_params(axis="x", pad=label_shift)
    ax.set_xlabel(str(x_name))
    ax.set_title(_title(settings) or f"{color_name} over {x_name}")
    ax.legend(title=str(color_name), bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.tight_layout()
    return fig


def _run_volcano_plot(adata, settings, de_info):
    group = _first_rank_group(adata, de_info, settings.get("de_volcano_group"))
    df = _rank_genes_groups_df(adata, group=group)
    if "logfoldchanges" not in df.columns:
        raise ScanpyPlotError("Volcano plot requires logfoldchanges in the Scanpy DE result.")
    p_col = "pvals_adj" if "pvals_adj" in df.columns else "pvals"
    if p_col not in df.columns:
        raise ScanpyPlotError("Volcano plot requires p-values in the Scanpy DE result.")
    plot_df = df.copy()
    plot_df["logfoldchanges"] = pd.to_numeric(plot_df["logfoldchanges"], errors="coerce")
    plot_df[p_col] = pd.to_numeric(plot_df[p_col], errors="coerce")
    max_abs_logfc = _as_optional_float(settings.get("de_volcano_max_abs_logfc"))
    if max_abs_logfc is not None:
        plot_df = plot_df[plot_df["logfoldchanges"].abs() <= max_abs_logfc]
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges", p_col])
    if plot_df.empty:
        raise ScanpyPlotError("No finite DE values are available for the volcano plot.")
    min_positive = np.nextafter(0, 1)
    plot_df["minus_log10_p"] = -np.log10(np.clip(plot_df[p_col].to_numpy(dtype=float), min_positive, None))
    fdr = float(settings.get("de_volcano_fdr") or settings.get("de_volcano_padj") or 0.05)
    abs_logfc = float(settings.get("de_volcano_logfc") or 1.0)
    significant = (plot_df[p_col] <= fdr) & (np.abs(plot_df["logfoldchanges"]) >= abs_logfc)

    fig, ax = plt.subplots(figsize=_figsize(settings))
    ax.scatter(
        plot_df.loc[~significant, "logfoldchanges"],
        plot_df.loc[~significant, "minus_log10_p"],
        s=float(settings.get("de_volcano_point_size") or 8),
        alpha=float(settings.get("alpha") or 0.45),
        c="lightgray",
        edgecolors="none",
        label="not significant",
    )
    up = significant & (plot_df["logfoldchanges"] > 0)
    down = significant & (plot_df["logfoldchanges"] < 0)
    ax.scatter(
        plot_df.loc[up, "logfoldchanges"],
        plot_df.loc[up, "minus_log10_p"],
        s=float(settings.get("de_volcano_point_size") or 8),
        alpha=float(settings.get("alpha") or 0.75),
        c="red",
        edgecolors="none",
        label="up",
    )
    ax.scatter(
        plot_df.loc[down, "logfoldchanges"],
        plot_df.loc[down, "minus_log10_p"],
        s=float(settings.get("de_volcano_point_size") or 8),
        alpha=float(settings.get("alpha") or 0.75),
        c="dodgerblue",
        edgecolors="none",
        label="down",
    )
    ax.axhline(-np.log10(fdr), linestyle="--", linewidth=1, color="green")
    ax.axvline(abs_logfc, linestyle="--", linewidth=1, color="black")
    ax.axvline(-abs_logfc, linestyle="--", linewidth=1, color="black")
    label_top_n = int(settings.get("de_volcano_label_top_n") or 15)
    if label_top_n > 0 and "names" in plot_df.columns:
        label_df = plot_df.sort_values(p_col, ascending=True).head(label_top_n)
        for _, row in label_df.iterrows():
            ax.text(
                row["logfoldchanges"],
                row["minus_log10_p"],
                str(row["names"]),
                fontsize=max(float(settings.get("font_size") or 10) - 2, 5),
            )
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel(f"-log10({p_col})")
    y_max = _as_optional_float(settings.get("de_volcano_ymax"))
    if y_max is not None:
        ax.set_ylim(top=y_max)
    ax.set_title(_title(settings) or f"Volcano plot: {group or 'DE'}")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    return fig

def run_standard_scanpy_plot(adata, settings):
    plot_type = settings.get("plot_type") or "embedding"
    layer = _prepare_x_for_plotting(adata, settings)
    use_raw = settings.get("use_raw")
    if use_raw == "auto" or use_raw == "":
        use_raw = None
    elif use_raw is not None:
        use_raw = _as_bool(use_raw)

    common = _json_kwargs(settings.get("advanced_kwargs"))
    if layer is not None:
        common.setdefault("layer", layer)
    if use_raw is not None:
        common.setdefault("use_raw", use_raw)

    with _figure_context(settings):
        if plot_type == "embedding":
            basis = _validate_basis(adata, settings.get("basis"))
            color = _listify_color(settings.get("color"))
            color = _validate_genes(adata, color, required=False) if color else None
            kwargs = {
                "basis": basis,
                "color": color,
                "legend_loc": settings.get("legend_loc") or "right margin",
                "legend_fontsize": _as_optional_int(settings.get("legend_fontsize")),
                "legend_fontoutline": _as_optional_int(settings.get("legend_fontoutline")),
                "size": _as_optional_float(settings.get("size")),
                "alpha": _as_optional_float(settings.get("alpha")),
                "components": settings.get("components") or None,
                "groups": _listify_color(settings.get("groups")),
                "palette": settings.get("palette") or None,
                "color_map": settings.get("color_map") or None,
                "sort_order": _as_bool(settings.get("sort_order"), True),
                "frameon": _as_bool(settings.get("frameon"), True),
                "add_outline": _as_bool(settings.get("add_outline"), False),
                "edges": _as_bool(settings.get("edges"), False),
                "arrows": _as_bool(settings.get("arrows"), False),
                "na_color": settings.get("na_color") or None,
                "title": _title(settings),
                "show": False,
            }
            kwargs.update(_common_plot_scale_kwargs(settings))
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.embedding(adata, **kwargs)

        elif plot_type in ("dotplot", "matrixplot", "heatmap", "stacked_violin", "tracksplot"):
            genes = _validate_genes(adata, settings.get("genes"), required=True)
            groupby = _validate_groupby(adata, settings.get("groupby"), required=True)
            _ensure_categorical_without_unused(adata, groupby)
            kwargs = {
                "var_names": genes,
                "groupby": groupby,
                **_common_grouped_plot_kwargs(settings),
            }
            if plot_type == "dotplot":
                kwargs["dot_min"] = _as_optional_float(settings.get("dot_min"))
                kwargs["dot_max"] = _as_optional_float(settings.get("dot_max"))
                kwargs["smallest_dot"] = _as_optional_float(settings.get("smallest_dot"))
                kwargs["expression_cutoff"] = _as_optional_float(settings.get("expression_cutoff"))
                kwargs["mean_only_expressed"] = _as_bool(settings.get("mean_only_expressed"), False)
                kwargs["dendrogram"] = _as_bool(settings.get("dendrogram"), False)
                kwargs["cmap"] = settings.get("color_map") or None
                kwargs["colorbar_title"] = settings.get("colorbar_title") or None
                kwargs["size_title"] = settings.get("size_title") or None
                plot_func = sc.pl.dotplot
            elif plot_type == "matrixplot":
                kwargs["cmap"] = settings.get("color_map") or None
                kwargs["dendrogram"] = _as_bool(settings.get("dendrogram"), False)
                kwargs["colorbar_title"] = settings.get("colorbar_title") or None
                plot_func = sc.pl.matrixplot
            elif plot_type == "heatmap":
                kwargs["cmap"] = settings.get("color_map") or None
                kwargs["dendrogram"] = _as_bool(settings.get("dendrogram"), False)
                kwargs["show_gene_labels"] = _as_optional_bool_select(settings.get("show_gene_labels"))
                plot_func = sc.pl.heatmap
            elif plot_type == "stacked_violin":
                kwargs["dendrogram"] = _as_bool(settings.get("dendrogram"), False)
                kwargs["cmap"] = settings.get("color_map") or None
                kwargs["stripplot"] = _as_bool(settings.get("stripplot"), True)
                kwargs["jitter"] = _as_optional_float(settings.get("jitter"))
                plot_func = sc.pl.stacked_violin
            else:
                kwargs["dendrogram"] = _as_bool(settings.get("dendrogram"), False)
                plot_func = sc.pl.tracksplot
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = plot_func(adata, **kwargs)

        elif plot_type == "violin":
            genes = _validate_genes(adata, settings.get("genes"), required=True)
            groupby = _validate_groupby(adata, settings.get("groupby"), required=False)
            if groupby:
                _ensure_categorical_without_unused(adata, groupby)
            kwargs = {
                "keys": genes,
                "groupby": groupby,
                "jitter": _as_optional_float(settings.get("jitter")),
                "rotation": _as_optional_float(settings.get("rotation")),
                "stripplot": _as_bool(settings.get("stripplot"), True),
                "multi_panel": _as_bool(settings.get("multi_panel"), False),
                "density_norm": settings.get("density_norm") or None,
                "order": _as_optional_categories(settings.get("categories_order")),
                "size": _as_optional_float(settings.get("size")),
                "show": False,
            }
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.violin(adata, **kwargs)

        elif plot_type == "embedding_density":
            ret = _run_embedding_density_plot(adata, settings, common)

        elif plot_type == "density_scatter":
            ret = _run_density_scatter_plot(adata, settings, layer=layer, use_raw=use_raw)

        elif plot_type == "stacked_barplot":
            ret = _run_stacked_barplot(adata, settings, layer=layer, use_raw=use_raw)

        elif plot_type == "highest_expr_genes":
            kwargs = {
                "n_top": int(settings.get("n_top") or 20),
                "show": False,
            }
            kwargs.update(common)
            ret = sc.pl.highest_expr_genes(adata, **kwargs)

        elif plot_type == "highly_variable_genes":
            if "highly_variable" not in adata.var.columns:
                raise ScanpyPlotError("This dataset does not contain adata.var['highly_variable'].")
            kwargs = {"show": False}
            kwargs.update(common)
            ret = sc.pl.highly_variable_genes(adata, **kwargs)

        elif plot_type == "pca_variance_ratio":
            if "pca" not in adata.uns:
                raise ScanpyPlotError("This dataset does not contain adata.uns['pca'].")
            kwargs = {
                "n_pcs": int(settings.get("n_pcs") or 30),
                "show": False,
            }
            kwargs.update(common)
            ret = sc.pl.pca_variance_ratio(adata, **kwargs)

        elif plot_type == "paga":
            if "paga" not in adata.uns:
                raise ScanpyPlotError("This dataset does not contain adata.uns['paga'].")
            kwargs = {
                "color": settings.get("color") or None,
                "threshold": _as_optional_float(settings.get("threshold")),
                "title": _title(settings),
                "show": False,
            }
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.paga(adata, **kwargs)

        elif plot_type == "dendrogram":
            groupby = _validate_groupby(adata, settings.get("groupby"), required=True)
            _ensure_categorical_without_unused(adata, groupby)
            dendro_key = f"dendrogram_{groupby}"
            if dendro_key not in adata.uns:
                sc.tl.dendrogram(adata, groupby=groupby)
            kwargs = {"groupby": groupby, "show": False}
            kwargs.update(common)
            ret = sc.pl.dendrogram(adata, **kwargs)

        else:
            raise ScanpyPlotError(f"Unsupported Scanpy plot type: {plot_type}")

        fig = _figure_from_return(ret)
        _set_figure_size(fig, settings)
        return _build_response(fig, settings, f"cellxgene_{plot_type}")


def _parse_de_use_raw(settings):
    if settings.get("_de_force_use_raw_false"):
        return False
    use_raw = settings.get("use_raw")
    if use_raw == "auto" or use_raw == "":
        return None
    if use_raw is not None:
        return _as_bool(use_raw)
    return None


def _prepare_de(adata, settings):
    settings = dict(settings or {})
    gene_filter = apply_de_gene_exclusion(adata, settings)
    layer = _prepare_x_for_plotting(adata, settings)
    method = settings.get("method") or "wilcoxon"
    n_genes = int(settings.get("n_genes") or 25)
    use_raw = _parse_de_use_raw(settings)
    corr_method = settings.get("corr_method") or "benjamini-hochberg"
    de_mode = settings.get("de_mode") or "selection"

    if de_mode == "obs_groups":
        groupby = _validate_groupby(adata, settings.get("de_groupby"), required=True)
        available = _ensure_categorical_without_unused(adata, groupby)
        if len(available) < 2:
            raise ScanpyPlotError("Differential expression by obs groups requires at least two groups in the current cell set.")

        raw_groups = str(settings.get("de_groups") or "all").strip()
        if raw_groups.lower() == "all":
            groups = "all"
            groups_to_plot = None
            csv_group = None
        else:
            groups_to_plot = parse_group_list(raw_groups)
            missing = [group for group in groups_to_plot if group not in available]
            if missing:
                raise ScanpyPlotError(f"Unknown DE group(s): {', '.join(missing)}")
            if not groups_to_plot:
                raise ScanpyPlotError("Enter one or more DE groups, or use 'all'.")
            groups = groups_to_plot
            csv_group = groups_to_plot[0] if len(groups_to_plot) == 1 else None

        reference = str(settings.get("de_reference") or "rest").strip()
        if reference != "rest" and reference not in available:
            raise ScanpyPlotError(f"Unknown DE reference group: {reference}")
        if groups_to_plot is None:
            plot_groupby_groups = list(available)
        else:
            plot_groupby_groups = [str(group) for group in groups_to_plot]
            if reference != "rest" and reference not in plot_groupby_groups:
                plot_groupby_groups.append(reference)
    else:
        if DE_GROUP_OBS not in adata.obs:
            raise ScanpyPlotError("Differential expression requires Selection 1 and Selection 2.")
        available = _ensure_categorical_without_unused(adata, DE_GROUP_OBS)
        group_counts = adata.obs[DE_GROUP_OBS].value_counts()
        if group_counts.get("Selection 1", 0) < 2 or group_counts.get("Selection 2", 0) < 2:
            raise ScanpyPlotError("Differential expression requires at least two cells in Selection 1 and Selection 2.")
        if "Selection 1" not in available or "Selection 2" not in available:
            raise ScanpyPlotError("Differential expression requires Selection 1 and Selection 2.")
        groupby = DE_GROUP_OBS
        groups = ["Selection 1"]
        groups_to_plot = ["Selection 1"]
        plot_groupby_groups = ["Selection 1", "Selection 2"]
        reference = "Selection 2"
        csv_group = "Selection 1"

    kwargs = {
        "groupby": groupby,
        "groups": groups,
        "reference": reference,
        "method": method,
        "n_genes": n_genes,
        "use_raw": use_raw,
        "pts": _as_bool(settings.get("pts"), True),
        "corr_method": corr_method,
        "rankby_abs": _as_bool(settings.get("rankby_abs"), False),
        "tie_correct": _as_bool(settings.get("tie_correct"), False),
    }
    if layer is not None:
        kwargs["layer"] = layer
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    sc.tl.rank_genes_groups(adata, **kwargs)

    return {
        "n_genes": n_genes,
        "groupby": groupby,
        "groups_to_plot": groups_to_plot,
        "plot_groupby_groups": plot_groupby_groups,
        "csv_group": csv_group,
        "gene_filter": gene_filter,
    }


def _requested_de_plot_n_genes(settings, computed_n_genes):
    requested = int(settings.get("n_genes") or computed_n_genes)
    if requested > computed_n_genes:
        raise ScanpyPlotError(
            f"This cached DE result was computed for {computed_n_genes} genes. "
            "Run DE again to plot more genes."
        )
    return requested


def _available_rank_groups(adata):
    try:
        names = adata.uns["rank_genes_groups"].get("names")
    except Exception:
        return []
    try:
        if getattr(names, "dtype", None) is not None and names.dtype.names:
            return [str(name) for name in names.dtype.names]
    except Exception:
        pass
    return []


def _requested_rank_groups_for_plot(adata, settings, de_info):
    requested = str(settings.get("de_plot_groups") or "").strip()
    if not requested:
        return de_info.get("groups_to_plot")
    if requested.lower() == "all":
        return None
    groups = parse_group_list(requested)
    available = _available_rank_groups(adata)
    if available:
        missing = [group for group in groups if group not in available]
        if missing:
            raise ScanpyPlotError(f"Unknown DE result group(s): {', '.join(missing[:10])}")
    return groups or de_info.get("groups_to_plot")


def _requested_de_plot_groupby_categories(adata, plot_groupby, settings, de_info):
    requested = str(settings.get("de_plot_groupby_groups") or "").strip()
    if not requested or requested.lower() in ("all", "*"):
        return None

    available = _ensure_categorical_without_unused(adata, plot_groupby)
    if requested.lower() in (
        "de",
        "de groups",
        "de comparison groups",
        "comparison groups",
        "comparison",
    ):
        groups = de_info.get("plot_groupby_groups") or de_info.get("groups_to_plot")
        if not groups:
            groups = _available_rank_groups(adata)
        groups = [str(group) for group in groups]
    else:
        groups = parse_group_list(requested)

    if not groups:
        return None

    missing = [group for group in groups if group not in available]
    if missing:
        preview = ", ".join(missing[:10])
        raise ScanpyPlotError(
            f"Unknown DE plot groupby category for {plot_groupby!r}: {preview}. "
            "Use category names from the DE plot groupby annotation, or leave this box blank."
        )
    return groups


def _de_plot_adata_for_grouped_plot(adata, settings, de_info, plot_type):
    if plot_type not in DE_GROUPED_PLOT_IDS:
        return adata, None

    plot_groupby = settings.get("de_plot_groupby") or de_info.get("groupby")
    if not plot_groupby:
        return adata, None
    _validate_groupby(adata, plot_groupby, required=True)

    selected_categories = _requested_de_plot_groupby_categories(
        adata,
        plot_groupby,
        settings,
        de_info,
    )
    if not selected_categories:
        _ensure_categorical_without_unused(adata, plot_groupby)
        return adata, None

    values = adata.obs[plot_groupby].astype(str).to_numpy()
    mask = np.isin(values, np.asarray(selected_categories, dtype=object))
    if not np.any(mask):
        raise ScanpyPlotError(
            "DE plot groupby category filter removed all cells. "
            "Choose categories present in the current plot data."
        )

    plot_adata = adata[mask, :].copy()
    remaining = _ensure_categorical_without_unused(plot_adata, plot_groupby)
    if len(remaining) < 1:
        raise ScanpyPlotError("DE plot groupby category filter removed all groups.")
    return plot_adata, selected_categories


def _rank_plot_common_kwargs(settings, de_info, adata=None):
    groups_to_plot = (
        _requested_rank_groups_for_plot(adata, settings, de_info)
        if adata is not None
        else de_info.get("groups_to_plot")
    )
    kwargs = {
        "n_genes": _requested_de_plot_n_genes(settings, int(de_info["n_genes"])),
        "groups": groups_to_plot,
        "show": False,
    }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    return kwargs


def _rank_grouped_plot_kwargs(settings, de_info, adata=None, plot_category_order=None):
    kwargs = _rank_plot_common_kwargs(settings, de_info, adata=adata)
    plot_groupby = settings.get("de_plot_groupby") or de_info["groupby"]
    if adata is not None and plot_groupby in adata.obs:
        _ensure_categorical_without_unused(adata, plot_groupby)
    categories_order = _as_optional_categories(settings.get("categories_order"))
    if categories_order is None and plot_category_order:
        categories_order = plot_category_order
    kwargs.update(
        {
            "groupby": plot_groupby,
            "values_to_plot": settings.get("de_values_to_plot") or None,
            "min_logfoldchange": _as_optional_float(settings.get("de_min_logfoldchange")),
            "cmap": settings.get("color_map") or None,
            "vmin": _as_optional_float(settings.get("vmin")),
            "vmax": _as_optional_float(settings.get("vmax")),
            "vcenter": _as_optional_float(settings.get("vcenter")),
            "dendrogram": _as_bool(settings.get("de_dendrogram"), False),
            "categories_order": categories_order,
        }
    )
    return {k: v for k, v in kwargs.items() if v is not None}


def render_differential_expression_plot(adata, settings, de_info):
    plot_type = settings.get("de_plot_type") or "rank_genes_groups"
    common = _json_kwargs(settings.get("advanced_kwargs"))
    csv_text, csv_preview = _csv_from_rank_genes_groups(adata, group=de_info.get("csv_group"))
    plot_adata, plot_category_order = _de_plot_adata_for_grouped_plot(
        adata,
        settings,
        de_info,
        plot_type,
    )
    plot_groupby = settings.get("de_plot_groupby") or de_info.get("groupby")
    if plot_groupby and plot_groupby in plot_adata.obs:
        _ensure_categorical_without_unused(plot_adata, plot_groupby)

    with _figure_context(settings):
        if plot_type == "rank_genes_groups":
            kwargs = _rank_plot_common_kwargs(settings, de_info, adata=adata)
            kwargs["sharey"] = _as_bool(settings.get("sharey"), False)
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups(adata, **kwargs)
        elif plot_type == "rank_genes_groups_dotplot":
            kwargs = _rank_grouped_plot_kwargs(
                settings,
                de_info,
                adata=plot_adata,
                plot_category_order=plot_category_order,
            )
            kwargs["standard_scale"] = settings.get("standard_scale") or None
            kwargs["swap_axes"] = _as_bool(settings.get("swap_axes"), False)
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups_dotplot(plot_adata, **kwargs)
        elif plot_type == "rank_genes_groups_matrixplot":
            kwargs = _rank_grouped_plot_kwargs(
                settings,
                de_info,
                adata=plot_adata,
                plot_category_order=plot_category_order,
            )
            kwargs["standard_scale"] = settings.get("standard_scale") or None
            kwargs["swap_axes"] = _as_bool(settings.get("swap_axes"), False)
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups_matrixplot(plot_adata, **kwargs)
        elif plot_type == "rank_genes_groups_stacked_violin":
            kwargs = _rank_grouped_plot_kwargs(
                settings,
                de_info,
                adata=plot_adata,
                plot_category_order=plot_category_order,
            )
            kwargs["standard_scale"] = settings.get("standard_scale") or None
            kwargs["swap_axes"] = _as_bool(settings.get("swap_axes"), False)
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups_stacked_violin(plot_adata, **kwargs)
        elif plot_type == "rank_genes_groups_tracksplot":
            kwargs = _rank_grouped_plot_kwargs(
                settings,
                de_info,
                adata=plot_adata,
                plot_category_order=plot_category_order,
            )
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups_tracksplot(plot_adata, **kwargs)
        elif plot_type == "rank_genes_groups_heatmap":
            kwargs = _rank_grouped_plot_kwargs(
                settings,
                de_info,
                adata=plot_adata,
                plot_category_order=plot_category_order,
            )
            kwargs["standard_scale"] = settings.get("standard_scale") or None
            kwargs["swap_axes"] = _as_bool(settings.get("swap_axes"), False)
            kwargs["show_gene_labels"] = _as_optional_bool_select(settings.get("show_gene_labels"))
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups_heatmap(plot_adata, **kwargs)
        elif plot_type == "rank_genes_groups_violin":
            kwargs = _rank_plot_common_kwargs(settings, de_info, adata=adata)
            kwargs.update(common)
            ret = sc.pl.rank_genes_groups_violin(adata, **kwargs)
        elif plot_type in ("volcano", "de_volcano"):
            ret = _run_volcano_plot(adata, settings, de_info)
        else:
            raise ScanpyPlotError(f"Unsupported differential-expression plot type: {plot_type}")

        fig = _figure_from_return(ret)
        _set_figure_size(fig, settings)
        response = _build_response(
            fig,
            settings,
            f"cellxgene_de_{plot_type}",
            csv_text=csv_text,
            csv_preview=csv_preview,
        )
        response["de_gene_filter"] = de_info.get("gene_filter")
        if "csv_preview" in response:
            response["de_table_preview"] = response["csv_preview"]
        return response


def run_differential_expression_plot(adata, settings):
    de_info = _prepare_de(adata, settings)
    return render_differential_expression_plot(adata, settings, de_info), de_info


def run_differential_expression_plot_from_precomputed(adata, settings, de_info):
    if "rank_genes_groups" not in adata.uns:
        raise ScanpyPlotError("Cached differential-expression result is incomplete. Run DE again.")
    return render_differential_expression_plot(adata, settings, de_info)

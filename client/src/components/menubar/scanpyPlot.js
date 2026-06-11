import React from "react";
import { connect } from "react-redux";
import { AnchorButton, Button, Callout, Checkbox, Classes, Code, Dialog, FormGroup, H4, HTMLSelect, HTMLTable, InputGroup, Tab, Tabs, TextArea, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";

const GENE_GROUP_TABS = new Set([
  "dotplot",
  "matrixplot",
  "heatmap",
  "violin",
  "stacked_violin",
  "tracksplot",
]);

const DE_GROUPED_PLOTS = new Set([
  "rank_genes_groups_dotplot",
  "rank_genes_groups_matrixplot",
  "rank_genes_groups_stacked_violin",
  "rank_genes_groups_tracksplot",
  "rank_genes_groups_heatmap",
]);

const SCANPY_DOCS = {
  embedding: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding.html",
    example: '{"ncols": 2, "wspace": 0.35}',
  },
  embedding_density: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding_density.html",
    example: '{"group": "all", "bg_dotsize": 40, "fg_dotsize": 120}',
  },
  density_scatter: {
    url: "https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.hexbin.html",
    example: '{"linewidths": 0.1}',
  },
  stacked_barplot: {
    url: "https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.bar.html",
    example: '{"edgecolor": "white", "linewidth": 0.2}',
  },
  dotplot: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.dotplot.html",
    example: '{"var_group_rotation": 45, "colorbar_title": "mean expression"}',
  },
  matrixplot: {
    url: "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.matrixplot.html",
    example: '{"colorbar_title": "mean expression", "var_group_rotation": 45}',
  },
  heatmap: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.heatmap.html",
    example: '{"show_gene_labels": true, "var_group_rotation": 45}',
  },
  violin: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.violin.html",
    example: '{"density_norm": "width", "ylabel": "expression"}',
  },
  stacked_violin: {
    url: "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.stacked_violin.html",
    example: '{"yticklabels": true, "row_palette": "tab20"}',
  },
  tracksplot: {
    url: "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.tracksplot.html",
    example:
      '{"var_group_labels": ["T", "B"], "var_group_positions": [[0, 2], [3, 5]]}',
  },
  highest_expr_genes: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.highest_expr_genes.html",
    example: '{"gene_symbols": "gene_symbols"}',
  },
  highly_variable_genes: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.highly_variable_genes.html",
    example: '{"log": true}',
  },
  pca_variance_ratio: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.pca_variance_ratio.html",
    example: '{"log": true}',
  },
  paga: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.paga.html",
    example: '{"layout": "fa", "edge_width_scale": 0.8}',
  },
  dendrogram: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.dendrogram.html",
    example: '{"orientation": "left"}',
  },
  rank_genes_groups: {
    url: "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.rank_genes_groups.html",
    example: '{"fontsize": 8}',
  },
  rank_genes_groups_dotplot: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_dotplot.html",
    example:
      '{"values_to_plot": "logfoldchanges", "cmap": "bwr", "vmin": -4, "vmax": 4}',
  },
  rank_genes_groups_matrixplot: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_matrixplot.html",
    example:
      '{"values_to_plot": "logfoldchanges", "cmap": "bwr", "vmin": -4, "vmax": 4}',
  },
  rank_genes_groups_stacked_violin: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_stacked_violin.html",
    example: '{"cmap": "viridis_r", "swap_axes": true}',
  },
  rank_genes_groups_tracksplot: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_tracksplot.html",
    example: '{"dendrogram": false}',
  },
  rank_genes_groups_heatmap: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_heatmap.html",
    example: '{"swap_axes": true, "show_gene_labels": false}',
  },
  rank_genes_groups_violin: {
    url: "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_violin.html",
    example: '{"strip": false}',
  },
  volcano: {
    url: "https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html",
    example: '{"alpha": 0.6}',
  },
};

function categoricalObsColumns(schema) {
  return (schema?.annotations?.obs?.columns || []).filter(
    (column) => column.type === "categorical" || column.categories
  );
}

function firstCategoricalObs(schema) {
  const columns = categoricalObsColumns(schema);
  return columns.length ? columns[0].name : "";
}

function defaultBasis(schema) {
  const layouts = schema?.layout?.obs || [];
  if (layouts.find((layout) => layout.name === "umap")) return "umap";
  return layouts.length ? layouts[0].name : "";
}

function dataUrl(mime, base64) {
  return `data:${mime};base64,${base64}`;
}

function downloadBase64(base64, mime, filename) {
  const raw = window.atob(base64);
  const bytes = new Uint8Array(raw.length);
  for (let i = 0; i < raw.length; i += 1) {
    bytes[i] = raw.charCodeAt(i);
  }
  const blob = new Blob([bytes], { type: mime });
  const url = window.URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  window.setTimeout(() => window.URL.revokeObjectURL(url), 0);
}

@connect((state) => ({
  annoMatrix: state.annoMatrix,
  differential: state.differential,
  scanpyPlot: state.scanpyPlot,
}))
class ScanpyPlot extends React.PureComponent {
  constructor(props) {
    super(props);
    const { annoMatrix } = props;
    const schema = annoMatrix?.schema;
    this.state = {
      isOpen: false,
      activeTabId: "embedding",
      zoom: 1,
      settings: {
        format: "png",
        dpi: "200",
        width: "7",
        height: "5",
        font_size: "10",
        title: "",
        transparent: false,
        useSelection1: false,
        useSelection2: false,
        basis: defaultBasis(schema),
        color: "",
        color_map: "viridis",
        legend_loc: "right margin",
        legend_fontsize: "10",
        legend_fontoutline: "",
        size: "",
        alpha: "",
        components: "",
        groups: "",
        palette: "",
        na_color: "lightgray",
        vmin: "",
        vmax: "",
        vcenter: "",
        sort_order: true,
        frameon: true,
        add_outline: false,
        edges: false,
        arrows: false,
        genes: "",
        groupby: firstCategoricalObs(schema),
        density_groupby: "",
        density_group: "all",
        density_key: "",
        density_bg_dotsize: "80",
        density_fg_dotsize: "180",
        ncols: "4",
        density_scatter_genes: "",
        density_scatter_cutoff: "0",
        density_scatter_bins: "50",
        stacked_bar_x: firstCategoricalObs(schema),
        stacked_bar_color: firstCategoricalObs(schema),
        stacked_bar_mode: "count",
        stacked_bar_x_order: "",
        stacked_bar_color_order: "",
        stacked_bar_gene_cutoff: "0",
        stacked_bar_label_rotation: "60",
        stacked_bar_x_font_size: "10",
        stacked_bar_label_shift: "0",
        dot_min: "",
        dot_max: "",
        smallest_dot: "0",
        expression_cutoff: "0",
        mean_only_expressed: false,
        standard_scale: "",
        swap_axes: false,
        dendrogram: false,
        categories_order: "",
        num_categories: "7",
        colorbar_title: "",
        size_title: "",
        show_gene_labels: "auto",
        stripplot: true,
        multi_panel: false,
        density_norm: "width",
        jitter: "0.4",
        rotation: "90",
        n_top: "20",
        n_pcs: "30",
        threshold: "",
        qc_plot_type: "highest_expr_genes",
        use_raw: "auto",
        layer: "",
        log1p: false,
        advanced_kwargs: "{}",
        de_mode: "selection",
        de_groupby: firstCategoricalObs(schema),
        de_groups: "all",
        de_reference: "rest",
        de_plot_groupby: "",
        de_plot_groupby_groups: "",
        de_plot_groups: "",
        de_exclude_genes: "",
        de_exclude_case_sensitive: false,
        de_plot_type: "rank_genes_groups",
        method: "wilcoxon",
        corr_method: "benjamini-hochberg",
        n_genes: "25",
        rankby_abs: false,
        tie_correct: false,
        pts: true,
        sharey: false,
        de_dendrogram: false,
        de_values_to_plot: "",
        de_min_logfoldchange: "",
        de_volcano_group: "",
        de_volcano_fdr: "0.05",
        de_volcano_logfc: "1",
        de_volcano_label_top_n: "15",
        de_volcano_point_size: "8",
        de_volcano_label_size: "8",
        de_volcano_max_abs_logfc: "",
        de_volcano_ymax: "",
      },
    };
  }

  componentDidUpdate(prevProps) {
    const { annoMatrix } = this.props;
    if (annoMatrix !== prevProps.annoMatrix) {
      this.updateDefaultsFromSchema();
    }
  }

  updateDefaultsFromSchema = () => {
    const { annoMatrix } = this.props;
    const schema = annoMatrix?.schema;
    this.setState((state) => {
      const { settings } = state;
      return {
        settings: {
          ...settings,
          basis: settings.basis || defaultBasis(schema),
          groupby: settings.groupby || firstCategoricalObs(schema),
          stacked_bar_x: settings.stacked_bar_x || firstCategoricalObs(schema),
          stacked_bar_color:
            settings.stacked_bar_color || firstCategoricalObs(schema),
          de_groupby: settings.de_groupby || firstCategoricalObs(schema),
        },
      };
    });
  };

  open = () => {
    this.updateDefaultsFromSchema();
    this.setState({ isOpen: true });
  };

  close = () => {
    this.setState({ isOpen: false });
  };

  setSetting = (name, value) => {
    this.setState((state) => ({
      settings: {
        ...state.settings,
        [name]: value,
      },
    }));
  };

  setActiveTab = (activeTabId) => {
    this.setState({ activeTabId });
  };

  plotTypeForActiveTab = () => {
    const { activeTabId, settings } = this.state;
    if (activeTabId === "qc_tools") return settings.qc_plot_type;
    if (activeTabId === "de") return settings.de_plot_type;
    return activeTabId;
  };

  activePlotHelp = () => {
    const plotType = this.plotTypeForActiveTab();
    return SCANPY_DOCS[plotType] || SCANPY_DOCS.embedding;
  };

  run = () => {
    const { dispatch } = this.props;
    const { activeTabId, settings } = this.state;
    const mode = activeTabId === "de" ? "de" : "plot";
    const plotSettings = {
      ...settings,
      plot_type:
        activeTabId === "qc_tools" ? settings.qc_plot_type : activeTabId,
    };
    dispatch(
      actions.runScanpyPlotAction({
        mode,
        useSelection1: settings.useSelection1,
        useSelection2: settings.useSelection2,
        settings: plotSettings,
      })
    );
  };

  regenerateDEPlotOnly = () => {
    const { dispatch, scanpyPlot } = this.props;
    const { settings } = this.state;
    const deResultId =
      scanpyPlot.lastDeResultId || scanpyPlot.result?.de_result_id;
    if (!deResultId) return;

    dispatch(
      actions.runScanpyPlotAction({
        mode: "de_plot_only",
        deResultId,
        settings,
      })
    );
  };

  zoomIn = () => {
    this.setState((state) => ({ zoom: Math.min(state.zoom * 1.25, 5) }));
  };

  zoomOut = () => {
    this.setState((state) => ({ zoom: Math.max(state.zoom / 1.25, 0.25) }));
  };

  resetZoom = () => {
    this.setState({ zoom: 1 });
  };

  renderSelectOptions = (options) =>
    options.map((option) => {
      if (typeof option === "string") {
        return (
          <option key={option} value={option}>
            {option || "None"}
          </option>
        );
      }
      return (
        <option key={option.value} value={option.value}>
          {option.label}
        </option>
      );
    });

  renderGeneralSettings() {
    const { differential, annoMatrix } = this.props;
    const { settings } = this.state;
    const currentViewCount = annoMatrix?.nObs || 0;
    const selection1Count = differential.celllist1?.length || 0;
    const selection2Count = differential.celllist2?.length || 0;

    return (
      <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 12 }}>
        <FormGroup label="Cells">
          <Checkbox
            checked={settings.useSelection1}
            disabled={selection1Count === 0}
            label={`Use Selection 1 (${selection1Count.toLocaleString()} cells)`}
            onChange={(e) => this.setSetting("useSelection1", e.target.checked)}
          />
          <Checkbox
            checked={settings.useSelection2}
            disabled={selection2Count === 0}
            label={`Use Selection 2 (${selection2Count.toLocaleString()} cells)`}
            onChange={(e) => this.setSetting("useSelection2", e.target.checked)}
          />
          <div style={{ fontSize: 12, color: globals.darkGrey }}>
            If neither selection is checked, the plot uses the current view (
            {currentViewCount.toLocaleString()} cells).
          </div>
        </FormGroup>
        <FormGroup label="Export settings">
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "1fr 1fr 1fr",
              gap: 8,
            }}
          >
            <FormGroup label="Format">
              <HTMLSelect
                value={settings.format}
                onChange={(e) => this.setSetting("format", e.target.value)}
              >
                {this.renderSelectOptions(["png", "svg", "pdf"])}
              </HTMLSelect>
            </FormGroup>
            <FormGroup label="DPI">
              <InputGroup
                value={settings.dpi}
                placeholder="200"
                onChange={(e) => this.setSetting("dpi", e.target.value)}
              />
            </FormGroup>
            <FormGroup label="Font size">
              <InputGroup
                value={settings.font_size}
                placeholder="10"
                onChange={(e) => this.setSetting("font_size", e.target.value)}
              />
            </FormGroup>
          </div>
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "1fr 1fr 1fr",
              gap: 8,
            }}
          >
            <FormGroup label="Width inches">
              <InputGroup
                value={settings.width}
                placeholder="7"
                onChange={(e) => this.setSetting("width", e.target.value)}
              />
            </FormGroup>
            <FormGroup label="Height inches">
              <InputGroup
                value={settings.height}
                placeholder="5"
                onChange={(e) => this.setSetting("height", e.target.value)}
              />
            </FormGroup>
            <FormGroup label="Background">
              <Checkbox
                checked={settings.transparent}
                label="Transparent"
                onChange={(e) =>
                  this.setSetting("transparent", e.target.checked)
                }
              />
            </FormGroup>
          </div>
          <FormGroup label="Title">
            <InputGroup
              value={settings.title}
              placeholder="Optional plot title"
              onChange={(e) => this.setSetting("title", e.target.value)}
            />
          </FormGroup>
        </FormGroup>
      </div>
    );
  }

  renderCommonDataSettings() {
    const { settings } = this.state;
    const { annoMatrix } = this.props;
    const layerNames = Object.keys(annoMatrix?.schema?.layers || {});
    return (
      <div
        style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 8 }}
      >
        <FormGroup label="Use raw">
          <HTMLSelect
            value={settings.use_raw}
            onChange={(e) => this.setSetting("use_raw", e.target.value)}
          >
            {this.renderSelectOptions([
              { value: "auto", label: "Auto" },
              { value: "true", label: "True" },
              { value: "false", label: "False" },
            ])}
          </HTMLSelect>
        </FormGroup>
        <FormGroup label="Layer">
          <HTMLSelect
            value={settings.layer}
            onChange={(e) => this.setSetting("layer", e.target.value)}
          >
            {this.renderSelectOptions([
              { value: "", label: "X" },
              ...layerNames.map((name) => ({ value: name, label: name })),
            ])}
          </HTMLSelect>
        </FormGroup>
        <FormGroup label="Transform">
          <Checkbox
            checked={settings.log1p}
            label="Apply log1p before plotting"
            onChange={(e) => this.setSetting("log1p", e.target.checked)}
          />
        </FormGroup>
      </div>
    );
  }

  renderAdvancedSettings() {
    const { settings } = this.state;
    const help = this.activePlotHelp();
    return (
      <FormGroup
        label={
          <span>
            Advanced Scanpy kwargs JSON{" "}
            <a href={help.url} rel="noopener noreferrer" target="_blank">
              docs
            </a>
          </span>
        }
        helperText={
          <span>
            Optional. Extra keyword arguments are passed to the selected
            scanpy.pl function. show, save, ax, and return_fig are ignored.
            Example: <Code>{help.example}</Code>
          </span>
        }
      >
        <TextArea
          fill
          growVertically
          value={settings.advanced_kwargs}
          onChange={(e) => this.setSetting("advanced_kwargs", e.target.value)}
        />
      </FormGroup>
    );
  }

  renderEmbeddingTab() {
    const { annoMatrix } = this.props;
    const { settings } = this.state;
    const layouts = annoMatrix?.schema?.layout?.obs || [];
    return (
      <div>
        <div
          style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 12 }}
        >
          <FormGroup label="Embedding basis">
            <HTMLSelect
              value={settings.basis}
              onChange={(e) => this.setSetting("basis", e.target.value)}
            >
              {this.renderSelectOptions(layouts.map((layout) => layout.name))}
            </HTMLSelect>
          </FormGroup>
          <FormGroup
            label="Color"
            helperText="Gene names or obs annotations. Separate multiple values by comma, space, or newline."
          >
            <InputGroup
              value={settings.color}
              placeholder="MS4A1, CD3D, leiden"
              onChange={(e) => this.setSetting("color", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Legend location">
            <HTMLSelect
              value={settings.legend_loc}
              onChange={(e) => this.setSetting("legend_loc", e.target.value)}
            >
              {this.renderSelectOptions([
                "right margin",
                "on data",
                "on data export",
                "best",
                "none",
              ])}
            </HTMLSelect>
          </FormGroup>
          <FormGroup label="Legend fontsize">
            <InputGroup
              value={settings.legend_fontsize}
              onChange={(e) =>
                this.setSetting("legend_fontsize", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Legend outline">
            <InputGroup
              value={settings.legend_fontoutline}
              placeholder="auto"
              onChange={(e) =>
                this.setSetting("legend_fontoutline", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Point size">
            <InputGroup
              value={settings.size}
              placeholder="auto"
              onChange={(e) => this.setSetting("size", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Color map">
            <InputGroup
              value={settings.color_map}
              onChange={(e) => this.setSetting("color_map", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="Alpha">
            <InputGroup
              value={settings.alpha}
              placeholder="auto"
              onChange={(e) => this.setSetting("alpha", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="Components">
            <InputGroup
              value={settings.components}
              placeholder="1,2"
              onChange={(e) => this.setSetting("components", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="Groups">
            <InputGroup
              value={settings.groups}
              placeholder="category1, category2"
              onChange={(e) => this.setSetting("groups", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="vmin">
            <InputGroup
              value={settings.vmin}
              placeholder="auto"
              onChange={(e) => this.setSetting("vmin", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vmax">
            <InputGroup
              value={settings.vmax}
              placeholder="auto"
              onChange={(e) => this.setSetting("vmax", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vcenter">
            <InputGroup
              value={settings.vcenter}
              placeholder="auto"
              onChange={(e) => this.setSetting("vcenter", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="NA color">
            <InputGroup
              value={settings.na_color}
              onChange={(e) => this.setSetting("na_color", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <Checkbox
            checked={settings.frameon}
            label="Draw frame"
            onChange={(e) => this.setSetting("frameon", e.target.checked)}
          />
          <Checkbox
            checked={settings.sort_order}
            label="Sort points by color"
            onChange={(e) => this.setSetting("sort_order", e.target.checked)}
          />
          <Checkbox
            checked={settings.add_outline}
            label="Add outline"
            onChange={(e) => this.setSetting("add_outline", e.target.checked)}
          />
          <Checkbox
            checked={settings.edges}
            label="Show graph edges"
            onChange={(e) => this.setSetting("edges", e.target.checked)}
          />
          <Checkbox
            checked={settings.arrows}
            label="Show arrows"
            onChange={(e) => this.setSetting("arrows", e.target.checked)}
          />
        </div>
        {this.renderCommonDataSettings()}
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderGroupBySelect(settingName = "groupby", options = {}) {
    const { annoMatrix } = this.props;
    const { settings } = this.state;
    const columns = categoricalObsColumns(annoMatrix?.schema);
    const columnNames = columns.map((column) => column.name);
    const withSelection = columnNames.includes("cellxgene_plot_selection")
      ? columnNames
      : ["cellxgene_plot_selection", ...columnNames];
    const finalOptions = options.emptyLabel
      ? [{ value: "", label: options.emptyLabel }, ...withSelection]
      : withSelection;
    return (
      <HTMLSelect
        value={settings[settingName]}
        onChange={(e) => this.setSetting(settingName, e.target.value)}
      >
        {this.renderSelectOptions(finalOptions)}
      </HTMLSelect>
    );
  }

  renderGeneByGroupTab(plotType) {
    const { settings } = this.state;
    const showDotSettings = plotType === "dotplot";
    const showViolinSettings = plotType === "violin";
    const showMatrixColorSettings = [
      "dotplot",
      "matrixplot",
      "heatmap",
      "stacked_violin",
    ].includes(plotType);
    const showGeneLabelSettings = plotType === "heatmap";
    return (
      <div>
        <div
          style={{ display: "grid", gridTemplateColumns: "2fr 1fr", gap: 12 }}
        >
          <FormGroup label="Genes">
            <TextArea
              fill
              growVertically
              value={settings.genes}
              placeholder="MS4A1, CD79A, CD3D"
              onChange={(e) => this.setSetting("genes", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="Group by categorical annotation">
            {this.renderGroupBySelect()}
            <div
              style={{ fontSize: 12, color: globals.darkGrey, marginTop: 8 }}
            >
              Use cellxgene_plot_selection to group by Selection 1 / Selection
              2.
            </div>
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          {showDotSettings ? (
            <FormGroup label="dot_min">
              <InputGroup
                value={settings.dot_min}
                placeholder="auto"
                onChange={(e) => this.setSetting("dot_min", e.target.value)}
              />
            </FormGroup>
          ) : null}
          {showDotSettings ? (
            <FormGroup label="dot_max">
              <InputGroup
                value={settings.dot_max}
                placeholder="auto"
                onChange={(e) => this.setSetting("dot_max", e.target.value)}
              />
            </FormGroup>
          ) : null}
          {showDotSettings ? (
            <FormGroup label="smallest_dot">
              <InputGroup
                value={settings.smallest_dot}
                placeholder="0"
                onChange={(e) =>
                  this.setSetting("smallest_dot", e.target.value)
                }
              />
            </FormGroup>
          ) : null}
          <FormGroup label="standard_scale">
            <HTMLSelect
              value={settings.standard_scale}
              onChange={(e) =>
                this.setSetting("standard_scale", e.target.value)
              }
            >
              {this.renderSelectOptions([
                { value: "", label: "None" },
                { value: "var", label: "var" },
                { value: "group", label: "group" },
              ])}
            </HTMLSelect>
          </FormGroup>
          <FormGroup label="Color map">
            <InputGroup
              value={settings.color_map}
              onChange={(e) => this.setSetting("color_map", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          {showDotSettings ? (
            <FormGroup label="expression_cutoff">
              <InputGroup
                value={settings.expression_cutoff}
                placeholder="0"
                onChange={(e) =>
                  this.setSetting("expression_cutoff", e.target.value)
                }
              />
            </FormGroup>
          ) : null}
          <FormGroup label="vmin">
            <InputGroup
              value={settings.vmin}
              placeholder="auto"
              onChange={(e) => this.setSetting("vmin", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vmax">
            <InputGroup
              value={settings.vmax}
              placeholder="auto"
              onChange={(e) => this.setSetting("vmax", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vcenter">
            <InputGroup
              value={settings.vcenter}
              placeholder="auto"
              onChange={(e) => this.setSetting("vcenter", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="num_categories">
            <InputGroup
              value={settings.num_categories}
              onChange={(e) =>
                this.setSetting("num_categories", e.target.value)
              }
            />
          </FormGroup>
        </div>
        {showMatrixColorSettings ? (
          <div
            style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 8 }}
          >
            <FormGroup label="Colorbar title">
              <InputGroup
                value={settings.colorbar_title}
                placeholder="auto"
                onChange={(e) =>
                  this.setSetting("colorbar_title", e.target.value)
                }
              />
            </FormGroup>
            {showDotSettings ? (
              <FormGroup label="Size title">
                <InputGroup
                  value={settings.size_title}
                  placeholder="auto"
                  onChange={(e) =>
                    this.setSetting("size_title", e.target.value)
                  }
                />
              </FormGroup>
            ) : null}
          </div>
        ) : null}
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <Checkbox
            checked={settings.swap_axes}
            label="Swap axes"
            onChange={(e) => this.setSetting("swap_axes", e.target.checked)}
          />
          <Checkbox
            checked={settings.dendrogram}
            label="Use / compute dendrogram when supported"
            onChange={(e) => this.setSetting("dendrogram", e.target.checked)}
          />
          {showDotSettings ? (
            <Checkbox
              checked={settings.mean_only_expressed}
              label="Mean only expressed"
              onChange={(e) =>
                this.setSetting("mean_only_expressed", e.target.checked)
              }
            />
          ) : null}
          {showGeneLabelSettings ? (
            <FormGroup label="Gene labels">
              <HTMLSelect
                value={settings.show_gene_labels}
                onChange={(e) =>
                  this.setSetting("show_gene_labels", e.target.value)
                }
              >
                {this.renderSelectOptions([
                  { value: "auto", label: "Auto" },
                  { value: "true", label: "True" },
                  { value: "false", label: "False" },
                ])}
              </HTMLSelect>
            </FormGroup>
          ) : null}
        </div>
        <FormGroup
          label="Categories order"
          helperText="Optional comma-separated order for the groupby categories."
        >
          <InputGroup
            value={settings.categories_order}
            placeholder="cluster 1, cluster 2, cluster 3"
            onChange={(e) =>
              this.setSetting("categories_order", e.target.value)
            }
          />
        </FormGroup>
        {showViolinSettings ? (
          <div>
            <div
              style={{
                display: "grid",
                gridTemplateColumns: "1fr 1fr 1fr",
                gap: 8,
              }}
            >
              <FormGroup label="Jitter">
                <InputGroup
                  value={settings.jitter}
                  onChange={(e) => this.setSetting("jitter", e.target.value)}
                />
              </FormGroup>
              <FormGroup label="Rotation">
                <InputGroup
                  value={settings.rotation}
                  onChange={(e) => this.setSetting("rotation", e.target.value)}
                />
              </FormGroup>
              <FormGroup label="Density norm">
                <HTMLSelect
                  value={settings.density_norm}
                  onChange={(e) =>
                    this.setSetting("density_norm", e.target.value)
                  }
                >
                  {this.renderSelectOptions(["width", "area", "count"])}
                </HTMLSelect>
              </FormGroup>
            </div>
            <Checkbox
              checked={settings.stripplot}
              label="Show stripplot"
              onChange={(e) => this.setSetting("stripplot", e.target.checked)}
            />
            <Checkbox
              checked={settings.multi_panel}
              label="Multi-panel"
              onChange={(e) => this.setSetting("multi_panel", e.target.checked)}
            />
          </div>
        ) : null}
        {this.renderCommonDataSettings()}
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderDensityTab() {
    const { annoMatrix } = this.props;
    const { settings } = this.state;
    const layouts = annoMatrix?.schema?.layout?.obs || [];
    return (
      <div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Embedding basis">
            <HTMLSelect
              value={settings.basis}
              onChange={(e) => this.setSetting("basis", e.target.value)}
            >
              {this.renderSelectOptions(layouts.map((layout) => layout.name))}
            </HTMLSelect>
          </FormGroup>
          <FormGroup
            label="Density groupby"
            helperText="Optional. If set, density is calculated per category."
          >
            {this.renderGroupBySelect("density_groupby", {
              emptyLabel: "Overall density",
            })}
          </FormGroup>
          <FormGroup label="Group to plot">
            <InputGroup
              value={settings.density_group}
              placeholder="all, none, or category names"
              onChange={(e) => this.setSetting("density_group", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Color map">
            <InputGroup
              value={settings.color_map}
              onChange={(e) => this.setSetting("color_map", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="Background dot size">
            <InputGroup
              value={settings.density_bg_dotsize}
              onChange={(e) =>
                this.setSetting("density_bg_dotsize", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Foreground dot size">
            <InputGroup
              value={settings.density_fg_dotsize}
              onChange={(e) =>
                this.setSetting("density_fg_dotsize", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Panels per row">
            <InputGroup
              value={settings.ncols}
              onChange={(e) => this.setSetting("ncols", e.target.value)}
            />
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="vmin">
            <InputGroup
              value={settings.vmin}
              placeholder="0"
              onChange={(e) => this.setSetting("vmin", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vmax">
            <InputGroup
              value={settings.vmax}
              placeholder="1"
              onChange={(e) => this.setSetting("vmax", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vcenter">
            <InputGroup
              value={settings.vcenter}
              placeholder="auto"
              onChange={(e) => this.setSetting("vcenter", e.target.value)}
            />
          </FormGroup>
        </div>
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderDensityScatterTab() {
    const { settings } = this.state;
    return (
      <div>
        <FormGroup
          label="Two genes or numeric annotations"
          helperText="Enter exactly two values. Example: CD83, CD33"
        >
          <InputGroup
            value={settings.density_scatter_genes}
            placeholder="CD83, CD33"
            onChange={(e) =>
              this.setSetting("density_scatter_genes", e.target.value)
            }
          />
        </FormGroup>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Expression cutoff line">
            <InputGroup
              value={settings.density_scatter_cutoff}
              onChange={(e) =>
                this.setSetting("density_scatter_cutoff", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Number of hexagons">
            <InputGroup
              value={settings.density_scatter_bins}
              onChange={(e) =>
                this.setSetting("density_scatter_bins", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Color map">
            <InputGroup
              value={settings.color_map}
              onChange={(e) => this.setSetting("color_map", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="Alpha">
            <InputGroup
              value={settings.alpha}
              placeholder="auto"
              onChange={(e) => this.setSetting("alpha", e.target.value)}
            />
          </FormGroup>
        </div>
        {this.renderCommonDataSettings()}
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderStackedBarplotTab() {
    const { settings } = this.state;
    return (
      <div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup
            label="X axis annotation or gene"
            helperText="Categorical obs annotation, or a gene to bin by expression cutoff."
          >
            <InputGroup
              value={settings.stacked_bar_x}
              placeholder="batch or CD83"
              onChange={(e) => this.setSetting("stacked_bar_x", e.target.value)}
            />
          </FormGroup>
          <FormGroup
            label="Color by annotation or gene"
            helperText="Categorical obs annotation, or a gene to bin by expression cutoff."
          >
            <InputGroup
              value={settings.stacked_bar_color}
              placeholder="cell_type or CD33"
              onChange={(e) =>
                this.setSetting("stacked_bar_color", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="Mode">
            <HTMLSelect
              value={settings.stacked_bar_mode}
              onChange={(e) =>
                this.setSetting("stacked_bar_mode", e.target.value)
              }
            >
              {this.renderSelectOptions([
                { value: "count", label: "Count" },
                { value: "proportion", label: "Proportion" },
                { value: "streamgraph", label: "Streamgraph" },
              ])}
            </HTMLSelect>
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Gene expression cutoff">
            <InputGroup
              value={settings.stacked_bar_gene_cutoff}
              onChange={(e) =>
                this.setSetting("stacked_bar_gene_cutoff", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="X label rotation">
            <InputGroup
              value={settings.stacked_bar_label_rotation}
              onChange={(e) =>
                this.setSetting("stacked_bar_label_rotation", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="X-axis font size">
            <InputGroup
              value={settings.stacked_bar_x_font_size}
              onChange={(e) =>
                this.setSetting("stacked_bar_x_font_size", e.target.value)
              }
            />
          </FormGroup>
          <FormGroup label="X-axis label shift">
            <InputGroup
              value={settings.stacked_bar_label_shift}
              onChange={(e) =>
                this.setSetting("stacked_bar_label_shift", e.target.value)
              }
            />
          </FormGroup>
        </div>
        <FormGroup
          label="X-axis category order"
          helperText="Optional comma-separated order. Categories not listed are appended."
        >
          <InputGroup
            value={settings.stacked_bar_x_order}
            placeholder="CAP1, CAP2, CAP3"
            onChange={(e) =>
              this.setSetting("stacked_bar_x_order", e.target.value)
            }
          />
        </FormGroup>
        <FormGroup
          label="Legend/color category order"
          helperText="Optional comma-separated order. Categories not listed are appended."
        >
          <InputGroup
            value={settings.stacked_bar_color_order}
            placeholder="B cells, T cells, NK cells"
            onChange={(e) =>
              this.setSetting("stacked_bar_color_order", e.target.value)
            }
          />
        </FormGroup>
        {this.renderCommonDataSettings()}
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderQCToolsTab() {
    const { settings } = this.state;
    return (
      <div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Plot">
            <HTMLSelect
              value={settings.qc_plot_type}
              onChange={(e) => this.setSetting("qc_plot_type", e.target.value)}
            >
              {this.renderSelectOptions([
                {
                  value: "highest_expr_genes",
                  label: "Highest expressed genes",
                },
                {
                  value: "highly_variable_genes",
                  label: "Highly variable genes",
                },
                { value: "pca_variance_ratio", label: "PCA variance ratio" },
                { value: "paga", label: "PAGA" },
                { value: "dendrogram", label: "Dendrogram" },
              ])}
            </HTMLSelect>
          </FormGroup>
          <FormGroup label="n_top / n_pcs">
            <InputGroup
              value={
                settings.qc_plot_type === "pca_variance_ratio"
                  ? settings.n_pcs
                  : settings.n_top
              }
              onChange={(e) => {
                const field =
                  settings.qc_plot_type === "pca_variance_ratio"
                    ? "n_pcs"
                    : "n_top";
                this.setSetting(field, e.target.value);
              }}
            />
          </FormGroup>
          <FormGroup label="Group by / color">
            {settings.qc_plot_type === "dendrogram" ? (
              this.renderGroupBySelect()
            ) : (
              <InputGroup
                value={settings.color}
                placeholder="annotation"
                onChange={(e) => this.setSetting("color", e.target.value)}
              />
            )}
          </FormGroup>
        </div>
        <FormGroup label="PAGA threshold">
          <InputGroup
            value={settings.threshold}
            placeholder="auto"
            onChange={(e) => this.setSetting("threshold", e.target.value)}
          />
        </FormGroup>
        {this.renderCommonDataSettings()}
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderDETab() {
    const { differential } = this.props;
    const { settings } = this.state;
    const selection1Count = differential.celllist1?.length || 0;
    const selection2Count = differential.celllist2?.length || 0;
    const selectionMode = settings.de_mode === "selection";
    const deGroupedPlot = DE_GROUPED_PLOTS.has(settings.de_plot_type);
    const volcanoPlot = settings.de_plot_type === "volcano";
    return (
      <div>
        <Callout
          intent={
            selectionMode && (!selection1Count || !selection2Count)
              ? "warning"
              : "primary"
          }
        >
          {selectionMode ? (
            <span>
              Differential expression compares Selection 1 (
              {selection1Count.toLocaleString()} cells) against Selection 2 (
              {selection2Count.toLocaleString()} cells). Set those selections
              with the existing 1 and 2 buttons before running this mode.
            </span>
          ) : (
            <span>
              Differential expression compares categories from the selected obs
              annotation. Use Groups = all for Scanpy multi-compare, or enter
              comma-separated groups.
            </span>
          )}
        </Callout>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr",
            gap: 8,
            marginTop: 12,
          }}
        >
          <FormGroup label="Comparison mode">
            <HTMLSelect
              value={settings.de_mode}
              onChange={(e) => this.setSetting("de_mode", e.target.value)}
            >
              {this.renderSelectOptions([
                { value: "selection", label: "Selection 1 vs Selection 2" },
                { value: "obs_groups", label: "Obs groups / multi-compare" },
              ])}
            </HTMLSelect>
          </FormGroup>
          <FormGroup label="Method">
            <HTMLSelect
              value={settings.method}
              onChange={(e) => this.setSetting("method", e.target.value)}
            >
              {this.renderSelectOptions([
                "wilcoxon",
                "t-test",
                "t-test_overestim_var",
                "logreg",
              ])}
            </HTMLSelect>
          </FormGroup>
          <FormGroup label="Correction">
            <HTMLSelect
              value={settings.corr_method}
              onChange={(e) => this.setSetting("corr_method", e.target.value)}
            >
              {this.renderSelectOptions(["benjamini-hochberg", "bonferroni"])}
            </HTMLSelect>
          </FormGroup>
        </div>
        {selectionMode ? null : (
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "1fr 1fr 1fr",
              gap: 8,
            }}
          >
            <FormGroup label="Obs groupby">
              {this.renderGroupBySelect("de_groupby")}
            </FormGroup>
            <FormGroup
              label="Groups"
              helperText="Use all, or comma-separated group names to test multiple groups."
            >
              <InputGroup
                value={settings.de_groups}
                placeholder="all or B cells, T cells"
                onChange={(e) => this.setSetting("de_groups", e.target.value)}
              />
            </FormGroup>
            <FormGroup label="Reference">
              <InputGroup
                value={settings.de_reference}
                placeholder="rest or one group name"
                onChange={(e) =>
                  this.setSetting("de_reference", e.target.value)
                }
              />
            </FormGroup>
          </div>
        )}
        <FormGroup
          label="Exclude genes from DE"
          helperText="Optional. Comma, space, or newline separated exact genes or wildcards. Example: MT-*, RPS*, RPL*. Matching is case-insensitive by default."
        >
          <TextArea
            fill
            growVertically
            value={settings.de_exclude_genes}
            placeholder="MT-*, RPS*, RPL*"
            onChange={(e) =>
              this.setSetting("de_exclude_genes", e.target.value)
            }
          />
          <Checkbox
            checked={settings.de_exclude_case_sensitive}
            label="Case-sensitive exclusion matching"
            onChange={(e) =>
              this.setSetting("de_exclude_case_sensitive", e.target.checked)
            }
          />
        </FormGroup>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="DE plot">
            <HTMLSelect
              value={settings.de_plot_type}
              onChange={(e) => this.setSetting("de_plot_type", e.target.value)}
            >
              {this.renderSelectOptions([
                { value: "rank_genes_groups", label: "Rank genes groups" },
                { value: "rank_genes_groups_dotplot", label: "Dot plot" },
                { value: "rank_genes_groups_matrixplot", label: "Matrix plot" },
                {
                  value: "rank_genes_groups_stacked_violin",
                  label: "Stacked violin",
                },
                { value: "rank_genes_groups_tracksplot", label: "Tracks plot" },
                { value: "rank_genes_groups_heatmap", label: "Heatmap" },
                { value: "rank_genes_groups_violin", label: "Violin" },
                { value: "volcano", label: "Volcano plot" },
              ])}
            </HTMLSelect>
          </FormGroup>
          <FormGroup label="Number of genes">
            <InputGroup
              value={settings.n_genes}
              onChange={(e) => this.setSetting("n_genes", e.target.value)}
            />
          </FormGroup>
          {deGroupedPlot ? (
            <FormGroup
              label="DE plot groupby"
              helperText="Controls grouping for the plot only; it does not change the DE comparison."
            >
              {this.renderGroupBySelect("de_plot_groupby", {
                emptyLabel: "Use DE comparison groups",
              })}
            </FormGroup>
          ) : null}
          {deGroupedPlot ? (
            <FormGroup
              label="Plot groupby categories to show"
              helperText="Optional. Blank/all = all categories; de = DE comparison groups; or enter comma-separated plot groupby categories."
            >
              <InputGroup
                value={settings.de_plot_groupby_groups}
                placeholder="all, de, or Treg, Th1"
                onChange={(e) =>
                  this.setSetting("de_plot_groupby_groups", e.target.value)
                }
              />
            </FormGroup>
          ) : null}
          <FormGroup label="Values to plot">
            <HTMLSelect
              value={settings.de_values_to_plot}
              onChange={(e) =>
                this.setSetting("de_values_to_plot", e.target.value)
              }
            >
              {this.renderSelectOptions([
                { value: "", label: "Expression / default" },
                { value: "scores", label: "scores" },
                { value: "logfoldchanges", label: "logfoldchanges" },
                { value: "pvals", label: "pvals" },
                { value: "pvals_adj", label: "pvals_adj" },
                { value: "log10_pvals", label: "log10_pvals" },
                { value: "log10_pvals_adj", label: "log10_pvals_adj" },
              ])}
            </HTMLSelect>
          </FormGroup>
        </div>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <FormGroup label="Color map">
            <InputGroup
              value={settings.color_map}
              onChange={(e) => this.setSetting("color_map", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vmin">
            <InputGroup
              value={settings.vmin}
              placeholder="auto"
              onChange={(e) => this.setSetting("vmin", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="vmax">
            <InputGroup
              value={settings.vmax}
              placeholder="auto"
              onChange={(e) => this.setSetting("vmax", e.target.value)}
            />
          </FormGroup>
          <FormGroup label="min_logfoldchange">
            <InputGroup
              value={settings.de_min_logfoldchange}
              placeholder="auto"
              onChange={(e) =>
                this.setSetting("de_min_logfoldchange", e.target.value)
              }
            />
          </FormGroup>
        </div>
        {volcanoPlot ? (
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "1fr 1fr 1fr 1fr",
              gap: 8,
            }}
          >
            <FormGroup label="Volcano group">
              <InputGroup
                value={settings.de_volcano_group}
                placeholder="auto or DE group name"
                onChange={(e) =>
                  this.setSetting("de_volcano_group", e.target.value)
                }
              />
            </FormGroup>
            <FormGroup label="Volcano FDR line">
              <InputGroup
                value={settings.de_volcano_fdr}
                onChange={(e) =>
                  this.setSetting("de_volcano_fdr", e.target.value)
                }
              />
            </FormGroup>
            <FormGroup label="Volcano abs logFC line">
              <InputGroup
                value={settings.de_volcano_logfc}
                onChange={(e) =>
                  this.setSetting("de_volcano_logfc", e.target.value)
                }
              />
            </FormGroup>
            <FormGroup label="Volcano labels / point size">
              <InputGroup
                value={`${settings.de_volcano_label_top_n}, ${settings.de_volcano_point_size}`}
                placeholder="15, 8"
                onChange={(e) => {
                  const [labels, size] = e.target.value.split(",");
                  this.setSetting("de_volcano_label_top_n", labels || "");
                  this.setSetting(
                    "de_volcano_point_size",
                    size === undefined ? "" : size.trim()
                  );
                }}
              />
            </FormGroup>
          </div>
        ) : null}
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "1fr 1fr 1fr 1fr",
            gap: 8,
          }}
        >
          <Checkbox
            checked={settings.sharey}
            label="Share y axis"
            onChange={(e) => this.setSetting("sharey", e.target.checked)}
          />
          <Checkbox
            checked={settings.rankby_abs}
            label="Rank by absolute score"
            onChange={(e) => this.setSetting("rankby_abs", e.target.checked)}
          />
          <Checkbox
            checked={settings.tie_correct}
            label="Tie correction"
            onChange={(e) => this.setSetting("tie_correct", e.target.checked)}
          />
          <Checkbox
            checked={settings.pts}
            label="Compute pct expressed"
            onChange={(e) => this.setSetting("pts", e.target.checked)}
          />
          <Checkbox
            checked={settings.de_dendrogram}
            label="Use dendrogram in DE plots"
            onChange={(e) => this.setSetting("de_dendrogram", e.target.checked)}
          />
          <Checkbox
            checked={settings.swap_axes}
            label="Swap axes"
            onChange={(e) => this.setSetting("swap_axes", e.target.checked)}
          />
        </div>
        <Callout intent="primary">
          Use <b>Run DE + plot</b> after changing DE method, comparison groups,
          reference, ranked gene count, raw/layer/log settings, or excluded
          genes. Use <b>Regenerate DE plot</b> when you only changed the DE plot
          type or plot-display options.
        </Callout>
        {this.renderCommonDataSettings()}
        {this.renderAdvancedSettings()}
      </div>
    );
  }

  renderTabPanel() {
    const { activeTabId } = this.state;
    if (activeTabId === "embedding") return this.renderEmbeddingTab();
    if (activeTabId === "embedding_density") return this.renderDensityTab();
    if (activeTabId === "density_scatter")
      return this.renderDensityScatterTab();
    if (activeTabId === "stacked_barplot")
      return this.renderStackedBarplotTab();
    if (GENE_GROUP_TABS.has(activeTabId))
      return this.renderGeneByGroupTab(activeTabId);
    if (activeTabId === "qc_tools") return this.renderQCToolsTab();
    if (activeTabId === "de") return this.renderDETab();
    return null;
  }

  renderCsvPreview() {
    const { scanpyPlot } = this.props;
    const preview = scanpyPlot.result?.csv_preview;
    if (!preview || !preview.columns || !preview.rows) return null;
    return (
      <div style={{ marginTop: 12 }}>
        <H4 style={{ marginBottom: 6 }}>DE CSV preview</H4>
        <div
          style={{
            maxHeight: 260,
            overflow: "auto",
            border: `1px solid ${globals.lightGrey}`,
            background: "white",
          }}
        >
          <HTMLTable className={Classes.HTML_TABLE} style={{ width: "100%" }}>
            <thead>
              <tr>
                {preview.columns.map((column) => (
                  <th key={column}>{column}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {preview.rows.map((row) => {
                const rowKey = preview.columns
                  .map((column) => String(row[column] || ""))
                  .join("|");
                return (
                  <tr key={rowKey}>
                    {preview.columns.map((column) => (
                      <td key={column}>
                        {row[column] === null || row[column] === undefined
                          ? ""
                          : String(row[column])}
                      </td>
                    ))}
                  </tr>
                );
              })}
            </tbody>
          </HTMLTable>
        </div>
        <div style={{ fontSize: 12, color: globals.darkGrey, marginTop: 4 }}>
          Showing {preview.rows.length.toLocaleString()} of{" "}
          {preview.row_count.toLocaleString()} rows.
        </div>
      </div>
    );
  }

  renderDETablePreview() {
    const { scanpyPlot } = this.props;
    const preview =
      scanpyPlot.result?.de_table_preview || scanpyPlot.result?.csv_preview;
    if (!preview) return null;

    let columns = [];
    let rows = [];
    let rowCount = 0;
    if (Array.isArray(preview)) {
      rows = preview;
      columns = rows.length ? Object.keys(rows[0]) : [];
      rowCount = rows.length;
    } else {
      columns = preview.columns || [];
      rows = preview.rows || [];
      rowCount = preview.row_count || rows.length;
    }

    if (!rows.length || !columns.length) return null;
    return (
      <div style={{ marginBottom: 12 }}>
        <b>DE result preview</b>{" "}
        <span style={{ color: globals.darkGrey }}>
          showing {rows.length.toLocaleString()} of {rowCount.toLocaleString()}{" "}
          rows
        </span>
        <div
          style={{
            maxHeight: 220,
            overflow: "auto",
            border: `1px solid ${globals.lightGrey}`,
            marginTop: 4,
          }}
        >
          <HTMLTable className={Classes.HTML_TABLE} style={{ width: "100%" }}>
            <thead>
              <tr>
                {columns.map((column) => (
                  <th key={column}>{column}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => {
                const rowKey = columns
                  .map((column) =>
                    row[column] === null || row[column] === undefined
                      ? ""
                      : String(row[column])
                  )
                  .join("|");
                return (
                  <tr key={rowKey}>
                    {columns.map((column) => (
                      <td key={column}>{String(row[column])}</td>
                    ))}
                  </tr>
                );
              })}
            </tbody>
          </HTMLTable>
        </div>
      </div>
    );
  }

  renderPreview() {
    const { scanpyPlot } = this.props;
    const { zoom } = this.state;
    const { result, loading, error } = scanpyPlot;
    if (loading) {
      return <Callout intent="primary">Running Scanpy plot...</Callout>;
    }
    if (error) {
      return <Callout intent="danger">{error}</Callout>;
    }
    if (!result) {
      return <Callout>Run a plot to preview it here.</Callout>;
    }

    return (
      <div>
        <div style={{ marginBottom: 8 }}>
          <Button small icon="zoom-in" onClick={this.zoomIn}>
            Zoom In
          </Button>{" "}
          <Button small icon="zoom-out" onClick={this.zoomOut}>
            Zoom Out
          </Button>{" "}
          <Button small onClick={this.resetZoom}>
            Reset
          </Button>{" "}
          <Button
            small
            icon="download"
            onClick={() =>
              downloadBase64(
                result.download_base64,
                result.download_mime,
                result.download_filename
              )
            }
          >
            Download plot
          </Button>{" "}
          {result.csv_base64 ? (
            <Button
              small
              icon="th"
              onClick={() =>
                downloadBase64(
                  result.csv_base64,
                  result.csv_mime,
                  result.csv_filename
                )
              }
            >
              Download DE CSV
            </Button>
          ) : null}
        </div>
        {this.renderDETablePreview()}
        <div
          style={{
            border: `1px solid ${globals.lightGrey}`,
            height: 420,
            overflow: "auto",
            background: "white",
          }}
        >
          <img
            src={dataUrl(result.preview_mime, result.preview_base64)}
            alt="Scanpy plot preview"
            style={{
              transform: `scale(${zoom})`,
              transformOrigin: "top left",
              maxWidth: "100%",
            }}
          />
        </div>
        {this.renderCsvPreview()}
      </div>
    );
  }

  renderDialog() {
    const { isOpen, activeTabId } = this.state;
    const { scanpyPlot } = this.props;
    return (
      <Dialog
        isOpen={isOpen}
        onClose={this.close}
        title="Create Scanpy plot"
        style={{ width: "min(1180px, 95vw)" }}
      >
        <div className={Classes.DIALOG_BODY}>
          {this.renderGeneralSettings()}
          <Tabs
            id="scanpy-plot-tabs"
            selectedTabId={activeTabId}
            onChange={this.setActiveTab}
            renderActiveTabPanelOnly
          >
            <Tab id="embedding" title="Embedding" />
            <Tab id="embedding_density" title="Density" />
            <Tab id="density_scatter" title="Density Scatter" />
            <Tab id="stacked_barplot" title="Stacked Barplot" />
            <Tab id="dotplot" title="DotPlot" />
            <Tab id="matrixplot" title="MatrixPlot" />
            <Tab id="heatmap" title="Heatmap" />
            <Tab id="violin" title="Violin" />
            <Tab id="stacked_violin" title="Stacked Violin" />
            <Tab id="tracksplot" title="TracksPlot" />
            <Tab id="qc_tools" title="QC / Tools" />
            <Tab id="de" title="Differential Expression" />
          </Tabs>
          <div style={{ marginTop: 12 }}>{this.renderTabPanel()}</div>
          <div style={{ marginTop: 12 }}>
            <Button
              intent="primary"
              icon="play"
              loading={scanpyPlot.loading}
              onClick={this.run}
            >
              {activeTabId === "de" ? "Run DE + plot" : "Run"}
            </Button>
            {activeTabId === "de" ? (
              <Button
                icon="refresh"
                disabled={scanpyPlot.loading || !scanpyPlot.lastDeResultId}
                style={{ marginLeft: 8 }}
                onClick={this.regenerateDEPlotOnly}
              >
                Regenerate DE plot
              </Button>
            ) : null}
          </div>
          <div style={{ marginTop: 16 }}>{this.renderPreview()}</div>
        </div>
      </Dialog>
    );
  }

  render() {
    return (
      <>
        <Tooltip
          content="Create Scanpy plots from the current view or saved selections"
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <AnchorButton
            className={styles.menubarButton}
            type="button"
            data-testid="scanpy-plot-button"
            icon="chart"
            onClick={this.open}
          />
        </Tooltip>
        {this.renderDialog()}
      </>
    );
  }
}

export default ScanpyPlot;

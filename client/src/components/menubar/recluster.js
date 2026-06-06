import React from "react";
import { connect } from "react-redux";
import { AnchorButton, Button, ButtonGroup, Callout, Checkbox, Classes, FormGroup, H5, HTMLSelect, NumericInput, Popover, Position, ProgressBar, TextArea, Tooltip } from "@blueprintjs/core";

import * as globals from "../../globals";
import styles from "./menubar.css";
import actions from "../../actions";

const DEFAULT_PARAMS = {
  resolution: 1.0,
  n_neighbors: 15,
  n_pcs: 50,
  use_rep: "X_pca",
  min_dist: 0.5,
  random_state: 0,
  gene_filter_mode: "none",
  gene_list: "",
  gene_filter_case_sensitive: false,
  gene_filter_log1p: false,
  gene_filter_scale: false,
};

function hasGeneTerms(text) {
  return !!text && text.split(/[\s,;]+/).some((term) => term.trim().length > 0);
}

@connect((state) => {
  const crossfilter = state.obsCrossfilter;
  const selectedCount = crossfilter?.countSelected?.() ?? 0;
  const limits = state.config?.limits ?? {};
  const representations =
    state.config?.parameters?.["recluster-representations"] ?? [];

  return {
    selectedCount,
    recluster: state.recluster,
    representations,
    minSelected: limits.recluster_cellcount_min ?? 10,
    maxSelected: limits.recluster_cellcount_max,
    maxGeneCount: limits.recluster_gene_count_max,
    maxExpressionValues: limits.recluster_expression_values_max,
    enabled: state.config?.parameters?.["recluster-enabled"] ?? false,
  };
})
class Recluster extends React.PureComponent {
  constructor(props) {
    super(props);
    const representations = props.representations.length
      ? props.representations
      : ["X_pca"];
    const useRep = representations.includes("X_pca")
      ? "X_pca"
      : representations[0];
    this.state = {
      params: {
        ...DEFAULT_PARAMS,
        use_rep: useRep,
      },
    };
  }

  componentDidUpdate(prevProps) {
    const { representations } = this.props;
    if (
      representations === prevProps.representations ||
      !representations.length
    ) {
      return;
    }

    const { params } = this.state;
    if (!representations.includes(params.use_rep)) {
      const useRep = representations.includes("X_pca")
        ? "X_pca"
        : representations[0];
      this.setState({ params: { ...params, use_rep: useRep } });
    }
  }

  setNumberParam = (name, value) => {
    if (!Number.isFinite(value)) return;
    this.setState((prev) => ({
      params: { ...prev.params, [name]: value },
    }));
  };

  setTextParam = (name, value) => {
    this.setState((prev) => ({
      params: { ...prev.params, [name]: value },
    }));
  };

  setBooleanParam = (name, value) => {
    this.setState((prev) => ({
      params: { ...prev.params, [name]: value },
    }));
  };

  setRepresentation = (event) => {
    const useRep = event.currentTarget.value;
    this.setState((prev) => ({
      params: { ...prev.params, use_rep: useRep },
    }));
  };

  setGeneFilterMode = (event) => {
    const geneFilterMode = event.currentTarget.value;
    this.setState((prev) => ({
      params: { ...prev.params, gene_filter_mode: geneFilterMode },
    }));
  };

  recluster = () => {
    const { dispatch } = this.props;
    const { params } = this.state;
    dispatch(actions.reclusterSelectionAction(params));
  };

  renderStatus() {
    const { recluster } = this.props;
    if (!recluster.loading) return null;

    return (
      <div style={{ marginBottom: 10 }}>
        <div style={{ marginBottom: 6 }}>
          {recluster.stage ?? "Reclustering..."}
          {recluster.nVars ? ` (${recluster.nVars} genes)` : ""}
        </div>
        <ProgressBar value={recluster.progress ?? 0} />
      </div>
    );
  }

  renderError() {
    const { recluster } = this.props;
    if (!recluster.error || recluster.loading) return null;

    const message =
      typeof recluster.error === "string"
        ? recluster.error
        : recluster.error.message || "Reclustering failed.";

    return (
      <Callout
        intent="danger"
        title="Reclustering failed"
        style={{ marginBottom: 10 }}
      >
        {message}
      </Callout>
    );
  }

  renderPopover() {
    const {
      selectedCount,
      recluster,
      representations,
      maxGeneCount,
      maxExpressionValues,
    } = this.props;
    const { params } = this.state;
    const repChoices = representations.length ? representations : ["X_pca"];
    const usingExpressionGenes = params.gene_filter_mode !== "none";
    const usingGeneList = ["whitelist", "blacklist"].includes(
      params.gene_filter_mode
    );
    const missingGeneList = usingGeneList && !hasGeneTerms(params.gene_list);

    return (
      <div style={{ width: 340, padding: 12 }}>
        <H5 className={Classes.HEADING}>Recluster selection</H5>
        <p style={{ marginTop: 0 }}>Selected cells: {selectedCount}</p>

        <FormGroup label="Representation">
          <HTMLSelect
            fill
            disabled={usingExpressionGenes}
            value={params.use_rep}
            onChange={this.setRepresentation}
          >
            {repChoices.map((rep) => (
              <option value={rep} key={rep}>
                {rep}
              </option>
            ))}
          </HTMLSelect>
        </FormGroup>

        <FormGroup label="Gene filter">
          <HTMLSelect
            fill
            value={params.gene_filter_mode}
            onChange={this.setGeneFilterMode}
          >
            <option value="none">Use representation</option>
            <option value="all">All genes: recompute PCA from X</option>
            <option value="whitelist">Whitelist: only listed genes</option>
            <option value="blacklist">Blacklist: exclude listed genes</option>
          </HTMLSelect>
        </FormGroup>

        {usingExpressionGenes ? (
          <>
            <Callout style={{ marginBottom: 10 }}>
              Expression-gene reclustering recomputes PCA from selected
              expression values in X. The representation above is ignored.
              All-gene mode uses every gene; whitelist and blacklist modes also
              support wildcards such as MT-*, RPL*, and RPS*.
            </Callout>

            {usingGeneList ? (
              <>
                <FormGroup
                  label={
                    params.gene_filter_mode === "whitelist"
                      ? "Whitelist genes"
                      : "Blacklist genes"
                  }
                  helperText="Paste gene symbols separated by newlines, spaces, commas, or semicolons."
                  intent={missingGeneList ? "danger" : undefined}
                >
                  <TextArea
                    fill
                    growVertically
                    rows={5}
                    value={params.gene_list}
                    placeholder={
                      params.gene_filter_mode === "blacklist"
                        ? "MT-*\nRPL*\nRPS*"
                        : "MS4A1\nCD79A\nCD74"
                    }
                    onChange={(event) =>
                      this.setTextParam("gene_list", event.currentTarget.value)
                    }
                  />
                </FormGroup>

                <Checkbox
                  checked={params.gene_filter_case_sensitive}
                  label="Case-sensitive gene matching"
                  onChange={(event) =>
                    this.setBooleanParam(
                      "gene_filter_case_sensitive",
                      event.currentTarget.checked
                    )
                  }
                />
              </>
            ) : null}
            <Checkbox
              checked={params.gene_filter_log1p}
              label="Apply log1p before PCA"
              onChange={(event) =>
                this.setBooleanParam(
                  "gene_filter_log1p",
                  event.currentTarget.checked
                )
              }
            />
            <Checkbox
              checked={params.gene_filter_scale}
              label="Scale genes before PCA; keeps sparse matrices sparse"
              onChange={(event) =>
                this.setBooleanParam(
                  "gene_filter_scale",
                  event.currentTarget.checked
                )
              }
            />

            <p style={{ marginTop: 4, marginBottom: 10, fontSize: 12 }}>
              Limits: {maxGeneCount ?? "unlimited"} genes;{" "}
              {maxExpressionValues ?? "unlimited"} cell x gene values.
            </p>
          </>
        ) : null}

        <FormGroup label="Resolution">
          <NumericInput
            fill
            min={0.01}
            stepSize={0.1}
            value={params.resolution}
            onValueChange={(value) => this.setNumberParam("resolution", value)}
          />
        </FormGroup>

        <FormGroup label="Neighbors">
          <NumericInput
            fill
            min={2}
            majorStepSize={10}
            stepSize={1}
            value={params.n_neighbors}
            onValueChange={(value) =>
              this.setNumberParam("n_neighbors", Math.round(value))
            }
          />
        </FormGroup>

        <FormGroup
          label={
            usingExpressionGenes
              ? "PCA components"
              : "PCs / representation dimensions"
          }
        >
          <NumericInput
            fill
            min={1}
            majorStepSize={10}
            stepSize={1}
            value={params.n_pcs}
            onValueChange={(value) =>
              this.setNumberParam("n_pcs", Math.round(value))
            }
          />
        </FormGroup>

        <FormGroup label="UMAP min_dist">
          <NumericInput
            fill
            min={0}
            stepSize={0.1}
            value={params.min_dist}
            onValueChange={(value) => this.setNumberParam("min_dist", value)}
          />
        </FormGroup>

        <FormGroup label="Random seed">
          <NumericInput
            fill
            min={0}
            stepSize={1}
            value={params.random_state}
            onValueChange={(value) =>
              this.setNumberParam("random_state", Math.round(value))
            }
          />
        </FormGroup>

        {this.renderStatus()}
        {this.renderError()}

        <Button
          intent="primary"
          loading={recluster.loading}
          disabled={recluster.loading || missingGeneList}
          onClick={this.recluster}
          text="Run reclustering"
        />
      </div>
    );
  }

  render() {
    const { selectedCount, minSelected, maxSelected, recluster, enabled } =
      this.props;

    const tooFew = selectedCount < minSelected;
    const tooMany = !!maxSelected && selectedCount > maxSelected;
    const disabled = !enabled || recluster.loading || tooFew || tooMany;

    let tooltip = "Subset and recluster selected cells";
    if (!enabled) tooltip = "Reclustering is not available for this dataset";
    else if (tooFew)
      tooltip = `Select at least ${minSelected} cells to recluster`;
    else if (tooMany) {
      tooltip = `Select no more than ${maxSelected} cells to recluster`;
    } else if (recluster.loading) {
      tooltip = recluster.stage ?? "Reclustering...";
    }

    return (
      <ButtonGroup className={styles.menubarButton}>
        <Popover
          content={this.renderPopover()}
          position={Position.BOTTOM}
          disabled={disabled && !recluster.loading}
        >
          <Tooltip
            content={tooltip}
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
          >
            <AnchorButton
              type="button"
              data-testid="recluster-button"
              disabled={disabled}
              loading={recluster.loading}
              icon="layout-hierarchy"
            />
          </Tooltip>
        </Popover>
      </ButtonGroup>
    );
  }
}

export default Recluster;

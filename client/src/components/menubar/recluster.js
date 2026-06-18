import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  Button,
  ButtonGroup,
  Callout,
  Checkbox,
  Classes,
  FormGroup,
  H5,
  HTMLSelect,
  NumericInput,
  Popover,
  Position,
  ProgressBar,
  TextArea,
  Tooltip,
} from "@blueprintjs/core";

import * as globals from "../../globals";
import styles from "./menubar.css";
import actions from "../../actions";

const DEFAULT_BLACKLIST_GENES = "MT-*, RPS*, RPL*";

const DEFAULT_PARAMS = {
  resolution: 1.0,
  n_neighbors: 15,
  n_pcs: 50,
  use_rep: "X_pca",
  min_dist: 0.5,
  random_state: 0,
  gene_filter_mode: "blacklist",
  gene_list: DEFAULT_BLACKLIST_GENES,
  gene_filter_case_sensitive: false,
  gene_filter_log1p: false,
  gene_filter_scale: false,
  harmony_enabled: false,
  harmony_batch_key: "",
  harmony_max_iter_harmony: 10,
  harmony_theta: "",
};

function hasGeneTerms(text) {
  return !!text && text.split(/[\s,;]+/).some((term) => term.trim().length > 0);
}

function isFiniteTextNumber(text) {
  if (text === null || text === undefined) return false;
  const trimmed = String(text).trim();
  if (!trimmed || trimmed === "." || trimmed === "-") return false;
  return Number.isFinite(Number(trimmed));
}

function textNumberValue(text) {
  return Number(String(text).trim());
}

@connect((state) => {
  const crossfilter = state.obsCrossfilter;
  const selectedCount = crossfilter?.countSelected?.() ?? 0;
  const limits = state.config?.limits ?? {};
  const representations =
    state.config?.parameters?.["recluster-representations"] ?? [];
  const batchKeys = state.config?.parameters?.["recluster-batch-keys"] ?? [];

  return {
    selectedCount,
    recluster: state.recluster,
    representations,
    batchKeys,
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
    const harmonyBatchKey = props.batchKeys.length ? props.batchKeys[0] : "";
    this.state = {
      params: {
        ...DEFAULT_PARAMS,
        use_rep: useRep,
        harmony_batch_key: harmonyBatchKey,
      },
      paramText: {
        resolution: String(DEFAULT_PARAMS.resolution),
        min_dist: String(DEFAULT_PARAMS.min_dist),
      },
    };
  }

  componentDidUpdate(prevProps) {
    const { representations, batchKeys } = this.props;
    const { params } = this.state;

    if (
      representations !== prevProps.representations &&
      representations.length &&
      !representations.includes(params.use_rep)
    ) {
      const useRep = representations.includes("X_pca")
        ? "X_pca"
        : representations[0];
      this.setState((prev) => ({
        params: { ...prev.params, use_rep: useRep },
      }));
    }

    if (
      batchKeys !== prevProps.batchKeys &&
      batchKeys.length &&
      !batchKeys.includes(params.harmony_batch_key)
    ) {
      this.setState((prev) => ({
        params: { ...prev.params, harmony_batch_key: batchKeys[0] },
      }));
    }
  }

  setNumberParam = (name, value) => {
    if (!Number.isFinite(value)) return;
    this.setState((prev) => ({
      params: { ...prev.params, [name]: value },
    }));
  };

  setDecimalParam = (name, value, valueAsString) => {
    const text = valueAsString === undefined ? String(value) : valueAsString;
    this.setState((prev) => {
      const nextParams = Number.isFinite(value)
        ? { ...prev.params, [name]: value }
        : prev.params;
      return {
        params: nextParams,
        paramText: { ...prev.paramText, [name]: text },
      };
    });
  };

  commitDecimalParam = (name, min) => {
    this.setState((prev) => {
      const text = prev.paramText[name];
      if (!isFiniteTextNumber(text)) {
        return {
          paramText: { ...prev.paramText, [name]: String(prev.params[name]) },
        };
      }

      const parsed = Math.max(min, textNumberValue(text));
      return {
        params: { ...prev.params, [name]: parsed },
        paramText: { ...prev.paramText, [name]: String(parsed) },
      };
    });
  };

  paramsForSubmit = () => {
    const { params, paramText } = this.state;
    const nextParams = { ...params };

    ["resolution", "min_dist"].forEach((name) => {
      if (isFiniteTextNumber(paramText[name])) {
        nextParams[name] = textNumberValue(paramText[name]);
      }
    });

    return nextParams;
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
    this.setState((prev) => {
      const nextParams = { ...prev.params, gene_filter_mode: geneFilterMode };
      if (
        geneFilterMode === "whitelist" &&
        prev.params.gene_filter_mode === "blacklist" &&
        prev.params.gene_list === DEFAULT_BLACKLIST_GENES
      ) {
        nextParams.gene_list = "";
      } else if (geneFilterMode === "blacklist" && !prev.params.gene_list) {
        nextParams.gene_list = DEFAULT_BLACKLIST_GENES;
      }
      return { params: nextParams };
    });
  };

  setHarmonyEnabled = (event) => {
    const enabled = event.currentTarget.checked;
    const { batchKeys } = this.props;
    this.setState((prev) => ({
      params: {
        ...prev.params,
        harmony_enabled: enabled,
        harmony_batch_key:
          enabled && !prev.params.harmony_batch_key && batchKeys.length
            ? batchKeys[0]
            : prev.params.harmony_batch_key,
      },
    }));
  };

  setHarmonyBatchKey = (event) => {
    const batchKey = event.currentTarget.value;
    this.setState((prev) => ({
      params: { ...prev.params, harmony_batch_key: batchKey },
    }));
  };

  recluster = () => {
    const { dispatch } = this.props;
    dispatch(actions.reclusterSelectionAction(this.paramsForSubmit()));
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
      batchKeys,
      maxGeneCount,
      maxExpressionValues,
    } = this.props;
    const { params, paramText } = this.state;
    const repChoices = representations.length ? representations : ["X_pca"];
    const usingExpressionGenes = params.gene_filter_mode !== "none";
    const usingGeneList = ["whitelist", "blacklist"].includes(
      params.gene_filter_mode
    );
    const missingGeneList = usingGeneList && !hasGeneTerms(params.gene_list);
    const invalidResolution =
      !isFiniteTextNumber(paramText.resolution) ||
      textNumberValue(paramText.resolution) < 0.01;
    const invalidMinDist =
      !isFiniteTextNumber(paramText.min_dist) ||
      textNumberValue(paramText.min_dist) < 0;
    const invalidHarmonyBatch =
      params.harmony_enabled && !params.harmony_batch_key;

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
                        ? DEFAULT_BLACKLIST_GENES
                        : "MS4A1, CD79A, CD74"
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

        <FormGroup label="Harmony batch correction">
          <Checkbox
            checked={params.harmony_enabled}
            label="Use Harmony before neighbors, Leiden, and UMAP"
            disabled={!batchKeys.length}
            onChange={this.setHarmonyEnabled}
          />
          {!batchKeys.length ? (
            <p style={{ marginTop: 4, fontSize: 12 }}>
              No categorical obs columns are available as Harmony batch keys.
            </p>
          ) : null}
        </FormGroup>

        {params.harmony_enabled ? (
          <>
            <FormGroup
              label="Harmony batch obs column"
              helperText="Harmony corrects the PCA/representation before the neighbor graph is computed."
              intent={invalidHarmonyBatch ? "danger" : undefined}
            >
              <HTMLSelect
                fill
                value={params.harmony_batch_key}
                onChange={this.setHarmonyBatchKey}
              >
                {batchKeys.map((key) => (
                  <option value={key} key={key}>
                    {key}
                  </option>
                ))}
              </HTMLSelect>
            </FormGroup>

            <FormGroup label="Harmony iterations">
              <NumericInput
                fill
                min={1}
                majorStepSize={10}
                stepSize={1}
                value={params.harmony_max_iter_harmony}
                onValueChange={(value) =>
                  this.setNumberParam(
                    "harmony_max_iter_harmony",
                    Math.round(value)
                  )
                }
              />
            </FormGroup>

            <FormGroup
              label="Harmony theta"
              helperText="Leave blank for the harmonypy default. Larger values encourage stronger batch mixing."
            >
              <input
                className={Classes.INPUT}
                style={{ width: "100%" }}
                value={params.harmony_theta}
                placeholder="default"
                onChange={(event) =>
                  this.setTextParam("harmony_theta", event.currentTarget.value)
                }
              />
            </FormGroup>
          </>
        ) : null}

        <FormGroup
          label="Resolution"
          intent={invalidResolution ? "danger" : undefined}
          helperText={
            invalidResolution ? "Enter a decimal number >= 0.01." : undefined
          }
        >
          <NumericInput
            fill
            allowNumericCharactersOnly={false}
            min={0.01}
            stepSize={0.1}
            value={paramText.resolution}
            onBlur={() => this.commitDecimalParam("resolution", 0.01)}
            onValueChange={(value, valueAsString) =>
              this.setDecimalParam("resolution", value, valueAsString)
            }
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

        <FormGroup
          label="UMAP min_dist"
          intent={invalidMinDist ? "danger" : undefined}
          helperText={
            invalidMinDist ? "Enter a decimal number >= 0." : undefined
          }
        >
          <NumericInput
            fill
            allowNumericCharactersOnly={false}
            min={0}
            stepSize={0.1}
            value={paramText.min_dist}
            onBlur={() => this.commitDecimalParam("min_dist", 0)}
            onValueChange={(value, valueAsString) =>
              this.setDecimalParam("min_dist", value, valueAsString)
            }
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
          disabled={
            recluster.loading ||
            missingGeneList ||
            invalidResolution ||
            invalidMinDist ||
            invalidHarmonyBatch
          }
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

import * as globals from "../globals";
import {
  doJsonRequest,
  dispatchNetworkErrorMessageToUser,
} from "../util/actionHelpers";
import { AnnoMatrixObsCrossfilter } from "../annoMatrix";
import { _switchEmbedding } from "./embedding";

const DEFAULT_BLACKLIST_GENES = "MT-*, RPS*, RPL*";

const DEFAULT_RECLUSTER_PARAMS = {
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

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

async function postReclusterJob(selectedObsIndices, params) {
  const url = `${globals.API.prefix}${globals.API.version}recluster/obs/jobs`;
  return doJsonRequest(url, {
    method: "POST",
    headers: new Headers({
      "Content-Type": "application/json",
    }),
    body: JSON.stringify({
      filter: { obs: { index: selectedObsIndices } },
      params,
    }),
  });
}

function pollReclusterJob(jobId, dispatch) {
  const url = `${globals.API.prefix}${globals.API.version}recluster/obs/jobs/${jobId}`;

  function pollOnce() {
    return doJsonRequest(url).then((job) => {
      dispatch({
        type: "recluster: status",
        jobId,
        status: job.status,
        stage: job.stage,
        progress: job.progress,
        nObs: job.n_obs,
        nVars: job.n_vars,
        geneFilter: job.gene_filter,
      });

      if (job.status === "complete") return job.result;
      if (job.status === "error") {
        throw new Error(job.error || "Reclustering failed.");
      }

      return sleep(750).then(pollOnce);
    });
  }

  return pollOnce();
}

export const reclusterSelectionAction =
  (params = {}) =>
  async (dispatch, getState) => {
    const state = getState();
    const { obsCrossfilter } = state;
    const selectedObsIndices = Array.from(obsCrossfilter.allSelectedLabels());
    const selectedCount = selectedObsIndices.length;
    const mergedParams = { ...DEFAULT_RECLUSTER_PARAMS, ...params };
    if (
      mergedParams.gene_filter_mode === "none" ||
      mergedParams.gene_filter_mode === "all"
    ) {
      mergedParams.gene_list = "";
    }

    if (selectedCount === 0) {
      dispatchNetworkErrorMessageToUser(
        "Select at least one cell before reclustering."
      );
      return;
    }

    dispatch({
      type: "recluster: started",
      selectedCount,
      stage: "Submitting recluster job",
      progress: 0,
    });

    try {
      const job = await postReclusterJob(selectedObsIndices, mergedParams);
      dispatch({
        type: "recluster: status",
        jobId: job.job_id,
        status: job.status,
        stage: job.stage,
        progress: job.progress,
        nObs: job.n_obs,
        nVars: job.n_vars,
        geneFilter: job.gene_filter,
      });

      const result = await pollReclusterJob(job.job_id, dispatch);

      const latestState = getState();
      const { obsCrossfilter: latestCrossfilter } = latestState;
      let { annoMatrix } = latestState;

      if (!annoMatrix.schema.layout.obsByName[result.layout.name]) {
        annoMatrix = annoMatrix.addEmbedding(result.layout);
      }
      if (!annoMatrix.schema.annotations.obsByName[result.cluster.name]) {
        annoMatrix = annoMatrix.addObsAnnotation(result.cluster);
      }

      const [nextAnnoMatrix] = await _switchEmbedding(
        annoMatrix,
        latestCrossfilter,
        result.layout.name
      );
      const nextObsCrossfilter = await new AnnoMatrixObsCrossfilter(
        nextAnnoMatrix
      ).select("emb", result.layout.name, {
        mode: "all",
      });

      dispatch({
        type: "recluster: success",
        result,
        annoMatrix: nextAnnoMatrix,
        obsCrossfilter: nextObsCrossfilter,
        layoutChoice: result.layout.name,
        clusterName: result.cluster.name,
      });
    } catch (error) {
      const message = error.message || "Reclustering failed.";
      if (!error.userMessageDisplayed) {
        dispatchNetworkErrorMessageToUser(message);
      }
      dispatch({
        type: "recluster: error",
        error: message,
      });
    }
  };

export default reclusterSelectionAction;

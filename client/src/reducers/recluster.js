const initialState = {
  loading: false,
  jobId: null,
  status: null,
  stage: null,
  progress: null,
  error: null,
  selectedCount: 0,
  nObs: null,
  nVars: null,
  geneFilter: null,
  lastResultId: null,
};

const Recluster = (state = initialState, action) => {
  switch (action.type) {
    case "recluster: started":
      return {
        ...state,
        loading: true,
        jobId: null,
        status: "submitting",
        stage: action.stage,
        progress: action.progress ?? 0,
        selectedCount: action.selectedCount,
        error: null,
      };

    case "recluster: status":
      return {
        ...state,
        loading: action.status !== "complete" && action.status !== "error",
        jobId: action.jobId,
        status: action.status,
        stage: action.stage,
        progress: action.progress,
        nObs: action.nObs,
        nVars: action.nVars,
        geneFilter: action.geneFilter,
      };

    case "recluster: success":
      return {
        ...state,
        loading: false,
        status: "complete",
        stage: "Complete",
        progress: 1,
        error: null,
        lastResultId: action.result.result_id,
      };

    case "recluster: error":
      return {
        ...state,
        loading: false,
        status: "error",
        error: action.error,
        stage: "Error",
        progress: 1,
      };

    default:
      return state;
  }
};

export default Recluster;

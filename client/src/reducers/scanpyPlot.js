const initialState = {
  loading: false,
  error: null,
  result: null,
  nObs: null,
  lastDeResultId: null,
};

export default function scanpyPlot(state = initialState, action) {
  switch (action.type) {
    case "scanpy plot: started":
      return {
        ...state,
        loading: true,
        error: null,
        nObs: action.nObs,
      };

    case "scanpy plot: success":
      return {
        ...state,
        loading: false,
        error: null,
        result: action.result,
        nObs: action.nObs,
        lastDeResultId: action.result.de_result_id || state.lastDeResultId,
      };

    case "scanpy plot: error":
      return {
        ...state,
        loading: false,
        error: action.error,
      };

    default:
      return state;
  }
}

const initialState = {
  loading: false,
  error: null,
  nObs: null,
  filename: null,
};

export default function exportH5ad(state = initialState, action) {
  switch (action.type) {
    case "export h5ad: started":
      return {
        ...state,
        loading: true,
        error: null,
        nObs: action.nObs,
        filename: null,
      };

    case "export h5ad: success":
      return {
        ...state,
        loading: false,
        error: null,
        nObs: action.nObs,
        filename: action.filename,
      };

    case "export h5ad: error":
      return {
        ...state,
        loading: false,
        error: action.message || action.error?.message,
      };

    default:
      return state;
  }
}

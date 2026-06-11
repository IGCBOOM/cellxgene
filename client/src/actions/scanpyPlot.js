import * as globals from "../globals";
import {
  doJsonRequest,
  dispatchNetworkErrorMessageToUser,
} from "../util/actionHelpers";
import {
  getReclusterResultIds,
  getWritableObsAnnotationsForView,
  rangeEncodeIndicesExclusive,
} from "./exportH5ad";

function unique(values) {
  return Array.from(new Set(values));
}

function intersectWithCurrentView(values, currentViewSet) {
  return values.filter((value) => currentViewSet.has(value));
}

function buildObsFilter(indices) {
  return {
    obs: {
      index: rangeEncodeIndicesExclusive(indices),
    },
  };
}

function getPlotObsIndices(annoMatrix, differential, options) {
  const currentView = Array.from(annoMatrix.rowIndex.labels());
  const currentViewSet = new Set(currentView);
  const selection1 = differential.celllist1
    ? intersectWithCurrentView(
        Array.from(differential.celllist1),
        currentViewSet
      )
    : [];
  const selection2 = differential.celllist2
    ? intersectWithCurrentView(
        Array.from(differential.celllist2),
        currentViewSet
      )
    : [];
  const useSelection1 = Boolean(options.useSelection1);
  const useSelection2 = Boolean(options.useSelection2);
  const settings = options.settings || {};

  if (options.mode === "de" && settings.de_mode !== "obs_groups") {
    return {
      currentView,
      selection1,
      selection2,
      obsIndices: unique([...selection1, ...selection2]),
    };
  }

  if (!useSelection1 && !useSelection2) {
    return {
      currentView,
      selection1,
      selection2,
      obsIndices: currentView,
    };
  }

  const obsIndices = unique([
    ...(useSelection1 ? selection1 : []),
    ...(useSelection2 ? selection2 : []),
  ]);
  return {
    currentView,
    selection1,
    selection2,
    obsIndices,
  };
}

export const runScanpyPlotAction =
  (options = {}) =>
  async (dispatch, getState) => {
    const { annoMatrix, differential, scanpyPlot } = getState();
    if (!annoMatrix) return;

    const settings = options.settings || {};
    const plotOnlyDE = options.mode === "de_plot_only";

    if (plotOnlyDE) {
      const deResultId = options.deResultId || scanpyPlot.lastDeResultId;
      if (!deResultId) {
        dispatchNetworkErrorMessageToUser(
          "No cached differential-expression result is available. Run DE first."
        );
        return;
      }

      dispatch({
        type: "scanpy plot: started",
        nObs: scanpyPlot.nObs,
      });

      try {
        const url = `${globals.API.prefix}${globals.API.version}plots/scanpy/run`;
        const result = await doJsonRequest(url, {
          method: "POST",
          headers: new Headers({
            "Content-Type": "application/json",
          }),
          body: JSON.stringify({
            mode: "de_plot_only",
            de_result_id: deResultId,
            settings,
          }),
        });

        dispatch({
          type: "scanpy plot: success",
          result,
          nObs: result.n_obs || scanpyPlot.nObs,
        });
      } catch (error) {
        const message = error.message || "Scanpy plot failed.";
        if (!error.userMessageDisplayed) {
          dispatchNetworkErrorMessageToUser(message);
        }
        dispatch({
          type: "scanpy plot: error",
          error: message,
        });
      }
      return;
    }

    const { obsIndices, selection1, selection2 } = getPlotObsIndices(
      annoMatrix,
      differential,
      options
    );

    if (
      options.mode === "de" &&
      settings.de_mode !== "obs_groups" &&
      (!selection1.length || !selection2.length)
    ) {
      dispatchNetworkErrorMessageToUser(
        "Differential expression plots require Selection 1 and Selection 2."
      );
      return;
    }

    if (obsIndices.length === 0) {
      dispatchNetworkErrorMessageToUser(
        "No cells are available for this plot."
      );
      return;
    }

    dispatch({
      type: "scanpy plot: started",
      nObs: obsIndices.length,
    });

    try {
      const obsAnnotations = await getWritableObsAnnotationsForView(annoMatrix);
      const body = {
        mode: options.mode === "de" ? "de" : "plot",
        filter: buildObsFilter(obsIndices),
        selection1: buildObsFilter(selection1),
        selection2: buildObsFilter(selection2),
        obs_annotations: obsAnnotations,
        recluster_result_ids: getReclusterResultIds(annoMatrix),
        settings,
      };

      const url = `${globals.API.prefix}${globals.API.version}plots/scanpy/run`;
      const result = await doJsonRequest(url, {
        method: "POST",
        headers: new Headers({
          "Content-Type": "application/json",
        }),
        body: JSON.stringify(body),
      });

      dispatch({
        type: "scanpy plot: success",
        result,
        nObs: obsIndices.length,
      });
    } catch (error) {
      const message = error.message || "Scanpy plot failed.";
      if (!error.userMessageDisplayed) {
        dispatchNetworkErrorMessageToUser(message);
      }
      dispatch({
        type: "scanpy plot: error",
        error: message,
      });
    }
  };

export default runScanpyPlotAction;

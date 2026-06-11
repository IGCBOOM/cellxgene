import * as globals from "../globals";
import {
  doFetch,
  dispatchNetworkErrorMessageToUser,
} from "../util/actionHelpers";

export function rangeEncodeIndicesExclusive(indices, minRangeLength = 8) {
  const result = [];
  let i = 0;

  while (i < indices.length) {
    const begin = indices[i];
    let current = begin;
    i += 1;

    while (i < indices.length && indices[i] === current + 1) {
      current = indices[i];
      i += 1;
    }

    if (current - begin + 1 >= minRangeLength) {
      result.push([begin, current + 1]);
    } else {
      for (let value = begin; value <= current; value += 1) {
        result.push(value);
      }
    }
  }

  return result;
}

function getFilenameFromContentDisposition(contentDisposition) {
  if (!contentDisposition) return null;

  const utf8Match = contentDisposition.match(/filename\*=UTF-8''([^;]+)/i);
  if (utf8Match) {
    return decodeURIComponent(utf8Match[1].replace(/['"]/g, ""));
  }

  const match = contentDisposition.match(/filename=([^;]+)/i);
  if (!match) return null;
  return match[1].trim().replace(/^"|"$/g, "");
}

function saveBlob(blob, filename) {
  const url = window.URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  window.setTimeout(() => window.URL.revokeObjectURL(url), 0);
}

function getWritableObsAnnotationSchemas(annoMatrix) {
  return annoMatrix.schema.annotations.obs.columns.filter((schema) =>
    Boolean(schema.writable)
  );
}

export function getReclusterResultIds(annoMatrix) {
  const resultIds = new Set();

  annoMatrix.schema.annotations.obs.columns.forEach((schema) => {
    if (schema.reclusterResultId) resultIds.add(schema.reclusterResultId);
  });

  annoMatrix.schema.layout.obs.forEach((schema) => {
    if (schema.reclusterResultId) resultIds.add(schema.reclusterResultId);
  });

  return Array.from(resultIds);
}

export async function getWritableObsAnnotationsForView(annoMatrix) {
  const writableSchemas = getWritableObsAnnotationSchemas(annoMatrix);
  if (writableSchemas.length === 0) return [];

  const names = writableSchemas.map((schema) => schema.name);
  const df = await annoMatrix.fetch("obs", names);

  return writableSchemas.map((schema) => ({
    name: schema.name,
    type: schema.type,
    categories: schema.categories || [],
    values: Array.from(df.col(schema.name).asArray()),
  }));
}

export const downloadCurrentViewH5ADAction =
  () => async (dispatch, getState) => {
    const { annoMatrix } = getState();
    if (!annoMatrix) return;

    const obsIndices = Array.from(annoMatrix.rowIndex.labels());
    if (obsIndices.length === 0) {
      dispatchNetworkErrorMessageToUser(
        "Current view contains no cells to export."
      );
      return;
    }

    dispatch({
      type: "export h5ad: started",
      nObs: obsIndices.length,
    });

    try {
      const obsAnnotations = await getWritableObsAnnotationsForView(annoMatrix);
      const body = {
        filter: {
          obs: {
            index: rangeEncodeIndicesExclusive(obsIndices),
          },
        },
        obs_annotations: obsAnnotations,
        recluster_result_ids: getReclusterResultIds(annoMatrix),
      };

      const url = `${globals.API.prefix}${globals.API.version}export/h5ad`;
      const res = await doFetch(url, {
        method: "POST",
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify(body),
      });
      const blob = await res.blob();
      const filename =
        getFilenameFromContentDisposition(
          res.headers.get("Content-Disposition")
        ) || "cellxgene_current_view.h5ad";
      saveBlob(blob, filename);

      dispatch({
        type: "export h5ad: success",
        nObs: obsIndices.length,
        filename,
      });
    } catch (error) {
      if (!error.userMessageDisplayed) {
        dispatchNetworkErrorMessageToUser(error.message);
      }
      dispatch({
        type: "export h5ad: error",
        error,
        message: error.message,
      });
    }
  };

export default downloadCurrentViewH5ADAction;

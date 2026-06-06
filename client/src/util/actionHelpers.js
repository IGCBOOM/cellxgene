import sortBy from "lodash.sortby";
/* XXX: cough, cough, ... */
import { postNetworkErrorToast } from "../components/framework/toasters";

/*
dispatch an action error to the user.   Currently we use
async toasts.
*/
let networkErrorToastKey = null;
export const dispatchNetworkErrorMessageToUser = (message) => {
  if (!networkErrorToastKey) {
    networkErrorToastKey = postNetworkErrorToast(message);
  } else {
    postNetworkErrorToast(message, networkErrorToastKey);
  }
};

/*
Catch unexpected errors and make sure we don't lose them!
*/
export function catchErrorsWrap(fn, dispatchToUser = false) {
  return (dispatch, getState) => {
    fn(dispatch, getState).catch((error) => {
      console.error(error);
      if (dispatchToUser) {
        dispatchNetworkErrorMessageToUser(error.message);
      }
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

function getHeader(headers, name) {
  if (!headers) return null;
  if (typeof headers.get === "function") return headers.get(name);
  return headers[name] || headers[name.toLowerCase()] || null;
}

function getErrorMessageFromPayload(payload) {
  if (!payload) return null;
  if (typeof payload === "string") return payload;

  const { error, message, detail, title } = payload;
  if (typeof error === "string") return error;
  if (error && typeof error.message === "string") return error.message;
  if (typeof message === "string") return message;
  if (typeof detail === "string") return detail;
  if (typeof title === "string") return title;
  return null;
}

async function readErrorMessage(res) {
  const fallback = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
  const contentType = res.headers.get("Content-Type") || "";

  try {
    if (contentType.includes("application/json")) {
      const payload = await res.clone().json();
      return getErrorMessageFromPayload(payload) || fallback;
    }

    if (!contentType.includes("text/html")) {
      const text = await res.clone().text();
      const trimmed = text.trim();
      if (trimmed.length) {
        return trimmed.length > 500 ? `${trimmed.slice(0, 500)}...` : trimmed;
      }
    }
  } catch (error) {
    // Fall through to the generic message.
  }

  return fallback;
}

/**
 * Wrapper to perform async fetch with some modest error handling
 * and decoding.  Arguments are identical to standard fetch.
 */
export const doFetch = async (url, init = {}) => {
  let res;
  try {
    // add defaults to the fetch init param.
    init = {
      method: "get",
      credentials: "include",
      ...init,
    };
    res = await fetch(url, init);
  } catch (e) {
    // fetch() only rejects for transport/network failures, not HTTP 4xx/5xx.
    const msg = "Unexpected HTTP error";
    dispatchNetworkErrorMessageToUser(msg);
    e.userMessageDisplayed = true;
    throw e;
  }

  const acceptType = getHeader(init.headers, "Accept");
  const contentType = res.headers.get("Content-Type") || "";
  if (res.ok && (!acceptType || contentType.includes(acceptType))) {
    return res;
  }

  const msg = await readErrorMessage(res);
  dispatchNetworkErrorMessageToUser(msg);
  const error = new Error(msg);
  error.status = res.status;
  error.statusText = res.statusText;
  error.url = url;
  error.userMessageDisplayed = true;
  throw error;
};

/*
Wrapper to perform an async fetch and JSON decode response.
*/
export const doJsonRequest = async (url, init = {}) => {
  const headers = new Headers(init.headers || {});
  headers.set("Accept", "application/json");
  const res = await doFetch(url, {
    ...init,
    headers,
  });
  return res.json();
};

/*
Wrapper to perform an async fetch for binary data.
*/
export const doBinaryRequest = async (url, init = {}) => {
  const headers = new Headers(init.headers || {});
  headers.set("Accept", "application/octet-stream");
  const res = await doFetch(url, {
    ...init,
    headers,
  });
  return res.arrayBuffer();
};

/*
This function "packs" filter index lists into the more efficient
"range" form specified in the REST 0.2 spec.

Specifically, it turns an array of indices [0, 1, 2, 10, 11, 14, ...]
into a form that encodes runs of consecutive numbers as [min, max].
Array may not be sorted, but will only contain uniq values.

Parameters:
   indices - input array of numbers (index)
   minRangeLength - hint, min range length before it is encoded into range format.
   sorted - boolean hint indicating array is presorted, ascending order

So [1, 2, 3, 4, 10, 11, 14] -> [ [1, 4], [10, 11], 14]
*/
export const rangeEncodeIndices = (
  indices,
  minRangeLength = 3,
  sorted = false
) => {
  if (indices.length === 0) {
    return indices;
  }

  if (!sorted) {
    indices = sortBy(indices);
  }

  const result = new Array(indices.length);
  let resultTail = 0;

  let i = 0;
  while (i < indices.length) {
    const begin = indices[i];
    let current;
    do {
      current = indices[i];
      i += 1;
    } while (i < indices.length && indices[i] === current + 1);

    if (current - begin + 1 >= minRangeLength) {
      result[resultTail] = [begin, current];
      resultTail += 1;
    } else {
      for (let j = begin; j <= current; j += 1, resultTail += 1) {
        result[resultTail] = j;
      }
    }
  }

  result.length = resultTail;
  return result;
};

import { createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";

import cascadeReducers from "./cascade";
import undoable from "./undoable";
import config from "./config";
import annoMatrix from "./annoMatrix";
import obsCrossfilter from "./obsCrossfilter";
import categoricalSelection from "./categoricalSelection";
import continuousSelection from "./continuousSelection";
import graphSelection from "./graphSelection";
import colors from "./colors";
import differential from "./differential";
import recluster from "./recluster";
import exportH5ad from "./exportH5ad";
import scanpyPlot from "./scanpyPlot";
import layoutChoice from "./layoutChoice";
import controls from "./controls";
import annotations from "./annotations";
import genesets from "./genesets";
import genesetsUI from "./genesetsUI";
import autosave from "./autosave";
import centroidLabels from "./centroidLabels";
import pointDialation from "./pointDilation";
import { gcMiddleware as annoMatrixGC } from "../annoMatrix";

import undoableConfig from "./undoableConfig";

const Reducer = undoable(
  cascadeReducers([
    ["config", config],
    ["annoMatrix", annoMatrix],
    ["obsCrossfilter", obsCrossfilter],
    ["annotations", annotations],
    ["genesets", genesets],
    ["genesetsUI", genesetsUI],
    ["layoutChoice", layoutChoice],
    ["categoricalSelection", categoricalSelection],
    ["continuousSelection", continuousSelection],
    ["graphSelection", graphSelection],
    ["colors", colors],
    ["controls", controls],
    ["differential", differential],
    ["recluster", recluster],
    ["exportH5ad", exportH5ad],
    ["scanpyPlot", scanpyPlot],
    ["centroidLabels", centroidLabels],
    ["pointDilation", pointDialation],
    ["autosave", autosave],
  ]),
  [
    "annoMatrix",
    "obsCrossfilter",
    "categoricalSelection",
    "continuousSelection",
    "graphSelection",
    "colors",
    "controls",
    "differential",
    "layoutChoice",
    "centroidLabels",
    "genesets",
    "annotations",
  ],
  undoableConfig
);

const store = createStore(Reducer, applyMiddleware(thunk, annoMatrixGC));

export default store;

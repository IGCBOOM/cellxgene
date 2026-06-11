import React from "react";
import { connect } from "react-redux";

import * as globals from "../../globals";
import Logo from "../framework/logo";
import Truncate from "../util/truncate";
import InformationMenu from "./infoMenu";

const DATASET_TITLE_FONT_SIZE = 14;
const PROCESS_MEMORY_POLL_MS = 2000;

@connect((state) => {
  const { corpora_props: corporaProps } = state.config;
  const correctVersion =
    ["1.0.0", "1.1.0"].indexOf(corporaProps?.version?.corpora_schema_version) >
    -1;
  return {
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    libraryVersions: state.config?.library_versions,
    aboutLink: state.config?.links?.["about-dataset"],
    tosURL: state.config?.parameters?.about_legal_tos,
    privacyURL: state.config?.parameters?.about_legal_privacy,
    title: correctVersion ? corporaProps?.title : undefined,
  };
})
class LeftSideBar extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      memoryLabel: null,
      memoryTitle: null,
    };
    this.memoryPollHandle = null;
    this.isComponentMounted = false;
  }

  componentDidMount() {
    this.isComponentMounted = true;
    this.fetchProcessMemory();
    this.memoryPollHandle = window.setInterval(
      this.fetchProcessMemory,
      PROCESS_MEMORY_POLL_MS
    );
  }

  componentWillUnmount() {
    this.isComponentMounted = false;
    if (this.memoryPollHandle) {
      window.clearInterval(this.memoryPollHandle);
      this.memoryPollHandle = null;
    }
  }

  fetchProcessMemory = () => {
    const url = `${globals.API.prefix}${globals.API.version}server/memory`;

    window
      .fetch(url, {
        cache: "no-store",
        credentials: "same-origin",
        headers: { Accept: "application/json" },
      })
      .then((response) => {
        if (!response.ok) throw new Error(response.statusText);
        return response.json();
      })
      .then((payload) => {
        if (!this.isComponentMounted) return;

        const processRss = payload?.process_rss_human || "unknown";
        const pid = payload?.pid ? `PID ${payload.pid}` : "server process";
        const systemUsed = payload?.system_used_human || "unknown";
        const systemTotal = payload?.system_total_human || "unknown";
        const systemPercent =
          typeof payload?.system_percent === "number"
            ? `${payload.system_percent.toFixed(1)}%`
            : "unknown";

        this.setState({
          memoryLabel: `RAM: ${processRss}`,
          memoryTitle: `${pid}; system RAM: ${systemUsed} / ${systemTotal} (${systemPercent})`,
        });
      })
      .catch(() => {
        // Memory telemetry is useful but should never interrupt the UI.
      });
  };

  render() {
    const {
      datasetTitle,
      libraryVersions,
      aboutLink,
      privacyURL,
      tosURL,
      dispatch,
      title,
    } = this.props;
    const { memoryLabel, memoryTitle } = this.state;

    return (
      <div
        style={{
          paddingLeft: 8,
          paddingTop: 8,
          width: globals.leftSidebarWidth,
          zIndex: 1,
          borderBottom: `1px solid ${globals.lighterGrey}`,
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
        }}
      >
        <div style={{ flex: "0 0 auto" }}>
          <Logo size={28} />
          <span
            style={{
              fontSize: 24,
              position: "relative",
              top: -6,
              fontWeight: "bold",
              marginLeft: 5,
              color: globals.logoColor,
              userSelect: "none",
            }}
          >
            cell
            <span
              style={{
                position: "relative",
                top: 1,
                fontWeight: 300,
                fontSize: 24,
              }}
            >
              ×
            </span>
            gene
          </span>
        </div>
        {memoryLabel ? (
          <span
            title={memoryTitle}
            style={{
              color: globals.darkGrey,
              flex: "0 0 auto",
              fontSize: 11,
              marginLeft: 8,
              marginRight: 8,
              position: "relative",
              top: -2,
              userSelect: "none",
              whiteSpace: "nowrap",
            }}
          >
            {memoryLabel}
          </span>
        ) : null}
        <div style={{ marginRight: 5, height: "100%" }}>
          <span
            minimal
            style={{
              fontSize: DATASET_TITLE_FONT_SIZE,
              padding: "5px 10px",
            }}
          >
            <Truncate>
              <span style={{ maxWidth: 155 }} data-testid="header">
                {title ?? datasetTitle}
              </span>
            </Truncate>
          </span>
          <InformationMenu
            {...{
              libraryVersions,
              aboutLink,
              tosURL,
              privacyURL,
              dispatch,
            }}
          />
        </div>
      </div>
    );
  }
}

export default LeftSideBar;

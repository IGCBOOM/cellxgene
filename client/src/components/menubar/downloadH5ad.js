import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import styles from "./menubar.css";
import * as globals from "../../globals";

const DownloadH5AD = React.memo((props) => {
  const { loading, nObs, handleDownload } = props;
  const cellText = nObs === 1 ? "1 cell" : `${nObs.toLocaleString()} cells`;

  return (
    <Tooltip
      content={`Download the current view as an h5ad file (${cellText}). Includes user annotations and generated recluster annotations.`}
      position="bottom"
      hoverOpenDelay={globals.tooltipHoverOpenDelay}
    >
      <AnchorButton
        className={styles.menubarButton}
        type="button"
        data-testid="download-current-view-h5ad-button"
        disabled={loading || nObs === 0}
        loading={loading}
        icon="download"
        onClick={handleDownload}
      />
    </Tooltip>
  );
});

export default DownloadH5AD;

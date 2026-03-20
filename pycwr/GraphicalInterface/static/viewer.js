const root = document.querySelector(".workspace");
const MIN_CANVAS_ZOOM = 1.0;
const MAX_CANVAS_ZOOM = 3.2;
const PAN_EDGE_PADDING = 36;

const state = {
  token: root.dataset.apiToken || "",
  directory: root.dataset.defaultDirectory || "",
  catalog: { stations: [] },
  stationFilter: "",
  stationId: "",
  metadata: null,
  playing: false,
  timer: null,
  renderToken: 0,
  ui: {
    controlsOpen: false,
    sectionOpen: false,
  },
  section: {
    mode: false,
    start: null,
    end: null,
    objectUrl: "",
    renderToken: 0,
  },
  canvas: {
    image: null,
    fitScale: 1,
    zoom: 1,
    offsetX: 0,
    offsetY: 0,
    dragging: false,
    dragMoved: false,
    dragStartX: 0,
    dragStartY: 0,
    dragOffsetX: 0,
    dragOffsetY: 0,
  },
};

const els = {
  directory: document.getElementById("directory-input"),
  directoryHint: document.getElementById("directory-hint"),
  scanButton: document.getElementById("scan-button"),
  saveLink: document.getElementById("save-link"),
  controlsToggleButton: document.getElementById("controls-toggle-button"),
  controlsCloseButton: document.getElementById("controls-close-button"),
  controlDrawer: document.getElementById("control-drawer"),
  drawerScrim: document.getElementById("drawer-scrim"),
  sectionToggleButton: document.getElementById("section-toggle-button"),
  sectionCloseButton: document.getElementById("section-close-button"),
  stationSearch: document.getElementById("station-search"),
  stationList: document.getElementById("station-list"),
  fieldChipStrip: document.getElementById("field-chip-strip"),
  sweepSelect: document.getElementById("sweep-select"),
  mapToggle: document.getElementById("map-toggle"),
  metadata: document.getElementById("metadata-list"),
  currentStation: document.getElementById("current-station"),
  currentTime: document.getElementById("current-time"),
  currentField: document.getElementById("current-field"),
  currentSweep: document.getElementById("current-sweep"),
  title: document.getElementById("viewer-title"),
  subtitle: document.getElementById("viewer-subtitle"),
  prevButton: document.getElementById("prev-button"),
  playButton: document.getElementById("play-button"),
  nextButton: document.getElementById("next-button"),
  speedSelect: document.getElementById("speed-select"),
  timelineList: document.getElementById("timeline-list"),
  status: document.getElementById("status-line"),
  canvas: document.getElementById("plot-canvas"),
  imageStage: document.querySelector(".image-stage"),
  canvasHint: document.getElementById("canvas-hint"),
  sectionFloatingHud: document.getElementById("section-floating-hud"),
  sectionFloatingStatus: document.getElementById("section-floating-status"),
  sectionPickButton: document.getElementById("section-pick-button"),
  sectionOpenButton: document.getElementById("section-open-button"),
  sectionClearInlineButton: document.getElementById("section-clear-inline-button"),
  sectionPanel: document.getElementById("section-panel"),
  sectionModeToggle: document.getElementById("section-mode-toggle"),
  sectionStartX: document.getElementById("section-start-x"),
  sectionStartY: document.getElementById("section-start-y"),
  sectionEndX: document.getElementById("section-end-x"),
  sectionEndY: document.getElementById("section-end-y"),
  sectionPickAgainButton: document.getElementById("section-pick-again-button"),
  sectionRenderButton: document.getElementById("section-render-button"),
  sectionClearButton: document.getElementById("section-clear-button"),
  sectionSaveLink: document.getElementById("section-save-link"),
  sectionSummary: document.getElementById("section-summary"),
  sectionImage: document.getElementById("section-image"),
  sectionPlaceholder: document.getElementById("section-placeholder"),
};

function setStatus(message, isError = false) {
  els.status.textContent = message;
  els.status.classList.toggle("is-error", isError);
}

function clearOptions(select) {
  while (select.firstChild) {
    select.removeChild(select.firstChild);
  }
}

function appendOption(select, value, label) {
  const option = document.createElement("option");
  option.value = value;
  option.textContent = label;
  select.appendChild(option);
}

function currentStation() {
  return state.catalog.stations.find((item) => item.station_id === state.stationId) || null;
}

function currentFile() {
  const station = currentStation();
  if (!station || !station.files.length) {
    return null;
  }
  const index = Number(els.timelineList.dataset.activeIndex || station.files.length - 1);
  if (index < 0 || index >= station.files.length) {
    return station.files[station.files.length - 1];
  }
  return station.files[index];
}

function currentTimelineIndex() {
  const station = currentStation();
  if (!station || !station.files.length) {
    return -1;
  }
  const index = Number(els.timelineList.dataset.activeIndex || station.files.length - 1);
  return Math.max(0, Math.min(index, station.files.length - 1));
}

function currentSweep() {
  return Number(els.sweepSelect.value || 0);
}

function currentSweepProfile() {
  if (!state.metadata) {
    return null;
  }
  return (state.metadata.sweep_profiles || []).find((item) => item.sweep === currentSweep()) || null;
}

function activeField() {
  return els.fieldChipStrip.dataset.currentField || "dBZ";
}

function currentFieldLegend() {
  if (!state.metadata) {
    return null;
  }
  const legends = state.metadata.field_legends || {};
  return legends[activeField()] || null;
}

function displayFieldName(field) {
  if (field === "HCL") {
    return "HCL / 水凝物分类";
  }
  return field;
}

function formatDisplayTime(isoString) {
  if (!isoString) {
    return "Unknown";
  }
  return isoString.replace("T", " ");
}

function formatSite(site) {
  if (!site) {
    return "-";
  }
  return `${site.longitude.toFixed(3)}°, ${site.latitude.toFixed(3)}°`;
}

function roundKm(value) {
  return Math.round(value * 10) / 10;
}

function formatPoint(point) {
  if (!point) {
    return "-";
  }
  return `(${point.x.toFixed(1)}, ${point.y.toFixed(1)}) km`;
}

function sectionLengthKm() {
  if (!state.section.start || !state.section.end) {
    return null;
  }
  const dx = state.section.end.x - state.section.start.x;
  const dy = state.section.end.y - state.section.start.y;
  return Math.hypot(dx, dy);
}

function sectionReady() {
  return Boolean(state.section.start && state.section.end && currentFile() && state.metadata);
}

function shouldRefreshSection() {
  return state.ui.sectionOpen && sectionReady() && !state.playing;
}

async function fetchJson(url) {
  const response = await fetch(url, {
    headers: {
      "X-Pycwr-Token": state.token,
    },
  });
  const data = await response.json();
  if (!response.ok || data.ok === false) {
    throw new Error(data.error || `Request failed: ${response.status}`);
  }
  return data;
}

function renderRows(container, rows) {
  container.innerHTML = "";
  rows.forEach(([label, value]) => {
    const dt = document.createElement("dt");
    dt.textContent = label;
    const dd = document.createElement("dd");
    dd.textContent = value;
    container.appendChild(dt);
    container.appendChild(dd);
  });
}

function syncDrawerState() {
  const anyOpen = state.ui.controlsOpen || state.ui.sectionOpen;
  els.controlDrawer.classList.toggle("is-open", state.ui.controlsOpen);
  els.controlDrawer.setAttribute("aria-hidden", String(!state.ui.controlsOpen));
  els.sectionPanel.classList.toggle("is-open", state.ui.sectionOpen);
  els.sectionPanel.setAttribute("aria-hidden", String(!state.ui.sectionOpen));
  els.drawerScrim.hidden = !anyOpen;
  els.controlsToggleButton.classList.toggle("is-active", state.ui.controlsOpen);
  els.sectionToggleButton.classList.toggle("is-active", state.ui.sectionOpen || state.section.mode || Boolean(state.section.objectUrl));
  els.sectionOpenButton.disabled = !Boolean(state.section.objectUrl);
}

function openControls(open) {
  state.ui.controlsOpen = open;
  if (open) {
    state.ui.sectionOpen = false;
  }
  syncDrawerState();
}

function openSectionPanel(open) {
  state.ui.sectionOpen = open;
  if (open) {
    state.ui.controlsOpen = false;
  }
  syncDrawerState();
}

function syncCanvasHint() {
  if (state.section.mode) {
    els.canvasHint.textContent = "Section mode: click start and end on Cartesian PPI";
    els.canvas.classList.add("is-section-mode");
    return;
  }
  els.canvasHint.textContent = "Wheel zoom | Drag pan | Double-click reset";
  els.canvas.classList.remove("is-section-mode");
}

function parseInputValue(input) {
  const raw = input.value.trim();
  if (!raw) {
    return null;
  }
  const value = Number(raw);
  return Number.isFinite(value) ? roundKm(value) : null;
}

function syncSectionInputs() {
  const start = state.section.start;
  const end = state.section.end;
  els.sectionStartX.value = start ? start.x.toFixed(1) : "";
  els.sectionStartY.value = start ? start.y.toFixed(1) : "";
  els.sectionEndX.value = end ? end.x.toFixed(1) : "";
  els.sectionEndY.value = end ? end.y.toFixed(1) : "";
}

function syncSectionSummary() {
  if (!state.section.start && !state.section.end) {
    els.sectionSummary.textContent = "Pick two points on the main PPI or enter Cartesian coordinates below.";
    els.sectionFloatingStatus.textContent = "Press Pick Line, then click the section start point on the radar image.";
    return;
  }
  if (state.section.start && !state.section.end) {
    els.sectionSummary.textContent = `Start ${formatPoint(state.section.start)} selected. Pick the section end point.`;
    els.sectionFloatingStatus.textContent = `Start ${formatPoint(state.section.start)} selected. Click the section end point.`;
    return;
  }
  const length = sectionLengthKm();
  els.sectionSummary.textContent = `Start ${formatPoint(state.section.start)} | End ${formatPoint(state.section.end)} | Length ${length.toFixed(1)} km`;
  els.sectionFloatingStatus.textContent = `Line ready: ${formatPoint(state.section.start)} -> ${formatPoint(state.section.end)} (${length.toFixed(1)} km).`;
}

function revokeSectionObjectUrl() {
  if (state.section.objectUrl) {
    URL.revokeObjectURL(state.section.objectUrl);
    state.section.objectUrl = "";
  }
}

function syncSectionPreview() {
  els.sectionPlaceholder.hidden = Boolean(state.section.objectUrl);
  els.sectionImage.hidden = !state.section.objectUrl;
}

function clearSectionImage() {
  revokeSectionObjectUrl();
  els.sectionImage.removeAttribute("src");
  els.sectionSaveLink.href = "#";
  els.sectionSaveLink.download = "pycwr-section.png";
  syncSectionPreview();
  syncDrawerState();
}

function setSectionMode(enabled) {
  state.section.mode = enabled;
  els.sectionModeToggle.checked = enabled;
  if (enabled && els.mapToggle.checked) {
    els.mapToggle.checked = false;
    refreshPlot();
    setStatus("Section picking uses Cartesian PPI view. Province overlay was turned off.");
  }
  syncCanvasHint();
  syncDrawerState();
  syncSectionSummary();
  drawCanvas();
}

function beginSectionPicking(reset = false) {
  if (reset) {
    state.section.start = null;
    state.section.end = null;
    clearSectionImage();
  }
  openControls(false);
  openSectionPanel(false);
  setSectionMode(true);
  syncSectionInputs();
  syncSectionSummary();
  setStatus("Section picking enabled. Click the start point.");
}

function renderStationList() {
  const previous = state.stationId;
  clearOptions(els.stationList);
  const filtered = state.catalog.stations.filter((station) => {
    if (!state.stationFilter) {
      return true;
    }
    const needle = state.stationFilter.toLowerCase();
    return station.station_id.toLowerCase().includes(needle) || station.station_name.toLowerCase().includes(needle);
  });
  filtered.forEach((station) => {
    appendOption(
      els.stationList,
      station.station_id,
      `${station.station_id} (${station.file_count})`
    );
  });
  if (!filtered.length) {
    appendOption(els.stationList, "", "No matching stations");
    els.stationList.disabled = true;
    return;
  }
  els.stationList.disabled = false;
  const nextStationId = filtered.some((item) => item.station_id === previous) ? previous : filtered[0].station_id;
  els.stationList.value = nextStationId;
  state.stationId = nextStationId;
}

function updateCurrentStrip() {
  const station = currentStation();
  const file = currentFile();
  const profile = currentSweepProfile();
  els.currentStation.textContent = station ? station.station_id : "-";
  els.currentTime.textContent = file ? formatDisplayTime(file.scan_time) : "-";
  els.currentField.textContent = displayFieldName(activeField());
  els.currentSweep.textContent = profile ? `${profile.fixed_angle.toFixed(1)} deg` : "-";
}

function renderFieldChips() {
  els.fieldChipStrip.innerHTML = "";
  if (!state.metadata || !state.metadata.fields) {
    return;
  }
  const currentField = activeField() || state.metadata.fields[0];
  els.fieldChipStrip.dataset.currentField = currentField;
  state.metadata.fields.forEach((field) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "field-chip";
    if (field === currentField) {
      button.classList.add("is-active");
    }
    button.textContent = displayFieldName(field);
    button.addEventListener("click", async () => {
      els.fieldChipStrip.dataset.currentField = field;
      renderFieldChips();
      updateCurrentStrip();
      updateTitle();
      refreshPlot();
      if (shouldRefreshSection()) {
        await refreshSection();
      }
    });
    els.fieldChipStrip.appendChild(button);
  });
}

function renderMetadata() {
  if (!state.metadata) {
    renderRows(els.metadata, []);
    return;
  }
  const rows = [
    ["Station", currentStation() ? currentStation().station_id : "-"],
    ["Scan", state.metadata.scan_type],
    ["Time", state.metadata.start_time],
    ["Site", formatSite(state.metadata.site)],
    ["Sweeps", String(state.metadata.sweeps.length)],
    ["Fields", state.metadata.fields.join(", ")],
  ];
  const legend = currentFieldLegend();
  if (legend && legend.length) {
    rows.push(["HCL Classes", legend.map((item) => `${item.value}:${item.label_zh}`).join(" | ")]);
  }
  renderRows(els.metadata, rows);
}

function renderTimeline() {
  els.timelineList.innerHTML = "";
  const station = currentStation();
  if (!station || !station.files.length) {
    return;
  }
  station.files.forEach((file, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "timeline-item";
    if (index === currentTimelineIndex()) {
      button.classList.add("is-active");
    }
    const time = file.scan_time ? file.scan_time.slice(11, 19) : "Unknown";
    const date = file.scan_time ? file.scan_time.slice(0, 10) : "";
    button.innerHTML = `<span class="timeline-item-time">${time}</span><span class="timeline-item-date">${date}</span>`;
    button.addEventListener("click", () => selectTimelineIndex(index));
    els.timelineList.appendChild(button);
  });
}

function updateTitle() {
  const station = currentStation();
  const file = currentFile();
  const profile = currentSweepProfile();
  if (!station || !file || !profile || !state.metadata) {
    els.title.textContent = "No file selected";
    els.subtitle.textContent = "Select a station and a timeline frame to load radar data.";
    return;
  }
  els.title.textContent = `${station.station_id} | ${activeField()} | Sweep ${profile.sweep}`;
  els.subtitle.textContent = `${formatDisplayTime(file.scan_time)} | ${state.metadata.scan_type.toUpperCase()} | ${formatSite(state.metadata.site)}`;
}

function plotUrl() {
  const file = currentFile();
  if (!file) {
    return "";
  }
  const stageRect = els.imageStage.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 1.6);
  const width = Math.max(760, Math.floor((stageRect.width - 12) * dpr));
  const height = Math.max(620, Math.floor((stageRect.height - 12) * dpr));
  const params = new URLSearchParams({
    token: state.token,
    path: file.path,
    field: activeField(),
    sweep: els.sweepSelect.value || "0",
    preset: "standard",
    continuous: "0",
    range_mode: "native",
    width: String(Math.min(width, 2200)),
    height: String(Math.min(height, 1800)),
  });
  const endpoint = els.mapToggle.checked ? "/plot/ppi_map.png" : "/plot/ppi.png";
  return `${endpoint}?${params.toString()}`;
}

function currentRangeKm() {
  const profile = currentSweepProfile();
  if (!profile) {
    return 150.0;
  }
  const fieldRanges = profile.field_ranges_km || {};
  const fieldRange = fieldRanges[activeField()];
  if (fieldRange && Number.isFinite(fieldRange.native) && fieldRange.native > 0) {
    return fieldRange.native;
  }
  if (Number.isFinite(profile.native_max_range_km) && profile.native_max_range_km > 0) {
    return profile.native_max_range_km;
  }
  if (Number.isFinite(profile.aligned_max_range_km) && profile.aligned_max_range_km > 0) {
    return profile.aligned_max_range_km;
  }
  return 150.0;
}

function currentPlotBounds() {
  const image = state.canvas.image;
  if (!image) {
    return null;
  }
  const leftPad = image.width * 0.07;
  const topPad = image.height * 0.065;
  const rightPad = image.width * 0.16;
  const bottomPad = image.height * 0.12;
  const usableWidth = Math.max(180, image.width - leftPad - rightPad);
  const usableHeight = Math.max(180, image.height - topPad - bottomPad);
  const side = Math.min(usableWidth, usableHeight);
  return {
    left: leftPad + Math.max(0, (usableWidth - side) / 2),
    top: topPad + Math.max(0, (usableHeight - side) / 2),
    width: side,
    height: side,
  };
}

function clamp(value, lower, upper) {
  return Math.min(Math.max(value, lower), upper);
}

function canvasMetrics() {
  const dpr = window.devicePixelRatio || 1;
  const width = Math.max(1, Math.floor(els.imageStage.clientWidth - 6));
  const height = Math.max(1, Math.floor(els.imageStage.clientHeight - 6));
  if (els.canvas.width !== Math.floor(width * dpr) || els.canvas.height !== Math.floor(height * dpr)) {
    els.canvas.width = Math.floor(width * dpr);
    els.canvas.height = Math.floor(height * dpr);
    els.canvas.style.width = `${width}px`;
    els.canvas.style.height = `${height}px`;
  }
  const ctx = els.canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  return { ctx, width, height };
}

function resetCanvasView() {
  const image = state.canvas.image;
  const bounds = currentPlotBounds();
  if (!image || !bounds) {
    return;
  }
  const { width, height } = canvasMetrics();
  const fitScale = 0.97 * Math.min(width / bounds.width, height / bounds.height);
  state.canvas.fitScale = fitScale;
  state.canvas.zoom = 1;
  state.canvas.offsetX = (width - bounds.width * fitScale) / 2 - bounds.left * fitScale;
  state.canvas.offsetY = (height - bounds.height * fitScale) / 2 - bounds.top * fitScale;
}

function constrainCanvasView() {
  const image = state.canvas.image;
  const bounds = currentPlotBounds();
  if (!image || !bounds) {
    return;
  }
  const width = Math.max(1, Math.floor(els.imageStage.clientWidth - 6));
  const height = Math.max(1, Math.floor(els.imageStage.clientHeight - 6));
  state.canvas.zoom = clamp(state.canvas.zoom, MIN_CANVAS_ZOOM, MAX_CANVAS_ZOOM);
  const scale = state.canvas.fitScale * state.canvas.zoom;
  const left = state.canvas.offsetX + bounds.left * scale;
  const right = state.canvas.offsetX + (bounds.left + bounds.width) * scale;
  const top = state.canvas.offsetY + bounds.top * scale;
  const bottom = state.canvas.offsetY + (bounds.top + bounds.height) * scale;
  const plotWidth = right - left;
  const plotHeight = bottom - top;

  if (plotWidth <= width - PAN_EDGE_PADDING * 2) {
    state.canvas.offsetX = (width - plotWidth) / 2 - bounds.left * scale;
  } else {
    const minOffsetX = width - PAN_EDGE_PADDING - (bounds.left + bounds.width) * scale;
    const maxOffsetX = PAN_EDGE_PADDING - bounds.left * scale;
    state.canvas.offsetX = clamp(state.canvas.offsetX, minOffsetX, maxOffsetX);
  }

  if (plotHeight <= height - PAN_EDGE_PADDING * 2) {
    state.canvas.offsetY = (height - plotHeight) / 2 - bounds.top * scale;
  } else {
    const minOffsetY = height - PAN_EDGE_PADDING - (bounds.top + bounds.height) * scale;
    const maxOffsetY = PAN_EDGE_PADDING - bounds.top * scale;
    state.canvas.offsetY = clamp(state.canvas.offsetY, minOffsetY, maxOffsetY);
  }
}

function radarKmToImage(point) {
  const bounds = currentPlotBounds();
  if (!bounds) {
    return null;
  }
  const rangeKm = currentRangeKm();
  if (!rangeKm || rangeKm <= 0) {
    return null;
  }
  return {
    x: bounds.left + ((point.x + rangeKm) / (2 * rangeKm)) * bounds.width,
    y: bounds.top + ((rangeKm - point.y) / (2 * rangeKm)) * bounds.height,
  };
}

function imageToRadarKm(point) {
  const bounds = currentPlotBounds();
  if (!bounds) {
    return null;
  }
  if (
    point.x < bounds.left ||
    point.x > bounds.left + bounds.width ||
    point.y < bounds.top ||
    point.y > bounds.top + bounds.height
  ) {
    return null;
  }
  const rangeKm = currentRangeKm();
  if (!rangeKm || rangeKm <= 0) {
    return null;
  }
  return {
    x: roundKm((((point.x - bounds.left) / bounds.width) * 2 - 1) * rangeKm),
    y: roundKm((1 - ((point.y - bounds.top) / bounds.height) * 2) * rangeKm),
  };
}

function imageToScreen(point) {
  const scale = state.canvas.fitScale * state.canvas.zoom;
  return {
    x: state.canvas.offsetX + point.x * scale,
    y: state.canvas.offsetY + point.y * scale,
  };
}

function screenToImage(x, y) {
  const scale = state.canvas.fitScale * state.canvas.zoom;
  return {
    x: (x - state.canvas.offsetX) / scale,
    y: (y - state.canvas.offsetY) / scale,
  };
}

function clientToRadarKm(clientX, clientY) {
  const rect = els.canvas.getBoundingClientRect();
  const imagePoint = screenToImage(clientX - rect.left, clientY - rect.top);
  return imageToRadarKm(imagePoint);
}

function drawSectionOverlay(ctx) {
  if (!state.canvas.image || (!state.section.mode && !state.section.start && !state.section.end)) {
    return;
  }
  const bounds = currentPlotBounds();
  if (!bounds) {
    return;
  }
  const topLeft = imageToScreen({ x: bounds.left, y: bounds.top });
  const bottomRight = imageToScreen({ x: bounds.left + bounds.width, y: bounds.top + bounds.height });
  ctx.save();
  ctx.strokeStyle = "rgba(38, 74, 112, 0.18)";
  ctx.lineWidth = 1;
  ctx.setLineDash([5, 5]);
  ctx.strokeRect(topLeft.x, topLeft.y, bottomRight.x - topLeft.x, bottomRight.y - topLeft.y);
  ctx.setLineDash([]);

  const startImage = state.section.start ? radarKmToImage(state.section.start) : null;
  const endImage = state.section.end ? radarKmToImage(state.section.end) : null;
  const start = startImage ? imageToScreen(startImage) : null;
  const end = endImage ? imageToScreen(endImage) : null;
  if (start && end) {
    ctx.strokeStyle = "#b33c2f";
    ctx.lineWidth = 2.3;
    ctx.beginPath();
    ctx.moveTo(start.x, start.y);
    ctx.lineTo(end.x, end.y);
    ctx.stroke();
  }
  [start, end].forEach((point, index) => {
    if (!point) {
      return;
    }
    ctx.fillStyle = index === 0 ? "#264a70" : "#b33c2f";
    ctx.strokeStyle = "#ffffff";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(point.x, point.y, 5.5, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
  });
  ctx.restore();
}

function drawCanvas() {
  const { ctx, width, height } = canvasMetrics();
  ctx.fillStyle = "#f4f7fb";
  ctx.fillRect(0, 0, width, height);
  const image = state.canvas.image;
  if (!image) {
    ctx.fillStyle = "#6b7785";
    ctx.font = "16px Times New Roman";
    ctx.fillText("Select station and timeline frame", 24, 32);
    return;
  }
  const scale = state.canvas.fitScale * state.canvas.zoom;
  ctx.save();
  ctx.translate(state.canvas.offsetX, state.canvas.offsetY);
  ctx.scale(scale, scale);
  ctx.drawImage(image, 0, 0);
  ctx.restore();
  drawSectionOverlay(ctx);
}

function applyZoom(clientX, clientY, factor) {
  if (!state.canvas.image) {
    return;
  }
  const bounds = currentPlotBounds();
  if (!bounds) {
    return;
  }
  const anchorImage = {
    x: bounds.left + bounds.width / 2,
    y: bounds.top + bounds.height / 2,
  };
  const anchorScreen = imageToScreen(anchorImage);
  state.canvas.zoom = clamp(state.canvas.zoom * factor, MIN_CANVAS_ZOOM, MAX_CANVAS_ZOOM);
  const scale = state.canvas.fitScale * state.canvas.zoom;
  state.canvas.offsetX = anchorScreen.x - anchorImage.x * scale;
  state.canvas.offsetY = anchorScreen.y - anchorImage.y * scale;
  constrainCanvasView();
  drawCanvas();
}

function setCanvasImage(url) {
  const token = ++state.renderToken;
  const image = new Image();
  image.onload = () => {
    if (token !== state.renderToken) {
      return;
    }
    state.canvas.image = image;
    resetCanvasView();
    drawCanvas();
  };
  image.onerror = () => {
    if (token !== state.renderToken) {
      return;
    }
    state.canvas.image = null;
    drawCanvas();
    setStatus("Unable to render radar image.", true);
  };
  image.src = `${url}&t=${Date.now()}`;
}

function refreshPlot() {
  const url = plotUrl();
  if (!url) {
    state.canvas.image = null;
    drawCanvas();
    return;
  }
  els.saveLink.href = url;
  els.saveLink.download = `${activeField()}_sweep${els.sweepSelect.value || 0}${els.mapToggle.checked ? "_map" : ""}.png`;
  updateCurrentStrip();
  updateTitle();
  setCanvasImage(url);
}

function sectionUrl() {
  const file = currentFile();
  if (!file || !sectionReady()) {
    return "";
  }
  const panelRect = els.sectionPanel.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 1.5);
  const width = Math.max(700, Math.floor(panelRect.width * dpr) || 1000);
  const height = Math.max(320, Math.floor((panelRect.height * 0.48) * dpr) || 460);
  const params = new URLSearchParams({
    token: state.token,
    path: file.path,
    field: activeField(),
    range_mode: "native",
    preset: "standard",
    width: String(Math.min(width, 2200)),
    height: String(Math.min(height, 1400)),
    start_x_km: state.section.start.x.toFixed(1),
    start_y_km: state.section.start.y.toFixed(1),
    end_x_km: state.section.end.x.toFixed(1),
    end_y_km: state.section.end.y.toFixed(1),
  });
  return `/plot/section.png?${params.toString()}`;
}

async function refreshSection() {
  if (!sectionReady()) {
    clearSectionImage();
    syncSectionSummary();
    return;
  }
  openSectionPanel(true);
  const url = sectionUrl();
  if (!url) {
    return;
  }
  const token = ++state.section.renderToken;
  els.sectionPlaceholder.hidden = false;
  els.sectionPlaceholder.textContent = "Rendering vertical section...";
  els.sectionImage.hidden = true;
  els.sectionSaveLink.href = url;
  els.sectionSaveLink.download = `${activeField()}_section.png`;
  try {
    const response = await fetch(url);
    const blob = await response.blob();
    if (token !== state.section.renderToken) {
      return;
    }
    revokeSectionObjectUrl();
    state.section.objectUrl = URL.createObjectURL(blob);
    els.sectionImage.src = state.section.objectUrl;
    syncSectionPreview();
    if (response.ok) {
      openSectionPanel(true);
      setStatus("Vertical section updated.");
    } else {
      openSectionPanel(true);
      setStatus("Section request returned an error preview.", true);
    }
  } catch (error) {
    if (token !== state.section.renderToken) {
      return;
    }
    clearSectionImage();
    els.sectionPlaceholder.hidden = false;
    els.sectionPlaceholder.textContent = "Unable to render section preview.";
    syncSectionPreview();
    setStatus(error.message, true);
  }
}

function readSectionInputs() {
  const startX = parseInputValue(els.sectionStartX);
  const startY = parseInputValue(els.sectionStartY);
  const endX = parseInputValue(els.sectionEndX);
  const endY = parseInputValue(els.sectionEndY);
  state.section.start = startX === null || startY === null ? null : { x: startX, y: startY };
  state.section.end = endX === null || endY === null ? null : { x: endX, y: endY };
  syncSectionInputs();
  syncSectionSummary();
  drawCanvas();
  if (!sectionReady()) {
    clearSectionImage();
  }
}

function resetSectionSelection() {
  state.section.start = null;
  state.section.end = null;
  setSectionMode(false);
  syncSectionInputs();
  syncSectionSummary();
  clearSectionImage();
  drawCanvas();
}

function pickSectionPoint(point) {
  if (!state.section.start || (state.section.start && state.section.end)) {
    state.section.start = point;
    state.section.end = null;
    clearSectionImage();
    setStatus("Section start selected. Click the end point.");
  } else {
    state.section.end = point;
    setSectionMode(false);
    setStatus("Rendering vertical section...");
    refreshSection();
  }
  syncSectionInputs();
  syncSectionSummary();
  drawCanvas();
}

async function loadMetadata(file) {
  if (!file) {
    return;
  }
  const payload = await fetchJson(`/api/metadata?path=${encodeURIComponent(file.path)}`);
  state.metadata = payload;
  const previousSweep = els.sweepSelect.value;
  clearOptions(els.sweepSelect);
  payload.sweeps.forEach((sweep) => appendOption(els.sweepSelect, String(sweep), String(sweep)));
  if (previousSweep && payload.sweeps.map(String).includes(previousSweep)) {
    els.sweepSelect.value = previousSweep;
  }
  if (!activeField() || !payload.fields.includes(activeField())) {
    els.fieldChipStrip.dataset.currentField = payload.fields[0];
  }
  renderFieldChips();
  renderMetadata();
  updateCurrentStrip();
  updateTitle();
  refreshPlot();
  syncSectionSummary();
  if (shouldRefreshSection()) {
    await refreshSection();
  }
}

async function selectTimelineIndex(index) {
  const station = currentStation();
  if (!station || index < 0 || index >= station.files.length) {
    return;
  }
  els.timelineList.dataset.activeIndex = String(index);
  renderTimeline();
  try {
    setStatus("Loading radar metadata...");
    await loadMetadata(station.files[index]);
    setStatus(`Loaded ${station.files[index].name}`);
  } catch (error) {
    setStatus(error.message, true);
  }
}

function stepTimeline(delta) {
  const station = currentStation();
  if (!station || !station.files.length) {
    return;
  }
  let nextIndex = currentTimelineIndex() + delta;
  if (nextIndex < 0) {
    nextIndex = station.files.length - 1;
  }
  if (nextIndex >= station.files.length) {
    nextIndex = 0;
  }
  selectTimelineIndex(nextIndex);
}

function togglePlayback() {
  if (state.playing) {
    window.clearInterval(state.timer);
    state.playing = false;
    els.playButton.textContent = "Play";
    setStatus("Playback stopped.");
    if (state.ui.sectionOpen && sectionReady()) {
      refreshSection();
    }
    return;
  }
  const speed = Number(els.speedSelect.value || "1");
  const intervalMs = Math.max(250, Math.floor(1000 / speed));
  state.playing = true;
  els.playButton.textContent = "Pause";
  setStatus("Playback started.");
  state.timer = window.setInterval(() => stepTimeline(1), intervalMs);
}

async function selectStation(stationId) {
  state.stationId = stationId;
  const station = currentStation();
  if (!station || !station.files.length) {
    renderTimeline();
    return;
  }
  els.timelineList.dataset.activeIndex = String(station.files.length - 1);
  renderTimeline();
  await selectTimelineIndex(station.files.length - 1);
}

function clearWorkspace() {
  state.catalog = { stations: [] };
  state.stationId = "";
  state.metadata = null;
  clearOptions(els.stationList);
  clearOptions(els.sweepSelect);
  els.fieldChipStrip.innerHTML = "";
  els.timelineList.innerHTML = "";
  renderRows(els.metadata, []);
  els.currentStation.textContent = "-";
  els.currentTime.textContent = "-";
  els.currentField.textContent = "-";
  els.currentSweep.textContent = "-";
  state.canvas.image = null;
  resetSectionSelection();
  drawCanvas();
  updateTitle();
}

async function scanDirectory() {
  try {
    setStatus("Building radar catalog...");
    const directory = els.directory.value.trim();
    const payload = await fetchJson(`/api/catalog?dir=${encodeURIComponent(directory)}`);
    state.directory = payload.directory;
    state.catalog = payload.catalog;
    els.directory.value = payload.directory;
    state.stationFilter = "";
    els.stationSearch.value = "";
    renderStationList();
    els.directoryHint.textContent = `${state.catalog.stations.length} stations discovered`;
    if (!state.catalog.stations.length) {
      clearWorkspace();
      setStatus("No supported radar files found.", true);
      return;
    }
    await selectStation(state.stationId);
  } catch (error) {
    clearWorkspace();
    setStatus(error.message, true);
  }
}

els.directory.value = state.directory;
els.playButton.textContent = "Play";
syncCanvasHint();
syncSectionInputs();
syncSectionSummary();
syncSectionPreview();
syncDrawerState();
drawCanvas();

els.controlsToggleButton.addEventListener("click", () => {
  openControls(!state.ui.controlsOpen);
});
els.controlsCloseButton.addEventListener("click", () => {
  openControls(false);
});
els.sectionToggleButton.addEventListener("click", () => {
  if (state.section.mode) {
    setSectionMode(false);
    setStatus("Section picking canceled.");
    return;
  }
  if (state.ui.sectionOpen) {
    openSectionPanel(false);
    return;
  }
  if (state.section.objectUrl) {
    openSectionPanel(true);
    return;
  }
  beginSectionPicking(true);
});
els.sectionCloseButton.addEventListener("click", () => {
  openSectionPanel(false);
});
els.drawerScrim.addEventListener("click", () => {
  openControls(false);
  openSectionPanel(false);
});

els.scanButton.addEventListener("click", scanDirectory);
els.stationSearch.addEventListener("input", (event) => {
  state.stationFilter = event.target.value.trim();
  renderStationList();
});
els.stationList.addEventListener("change", (event) => {
  if (!event.target.value) {
    return;
  }
  selectStation(event.target.value);
});
els.sweepSelect.addEventListener("change", async () => {
  refreshPlot();
  if (shouldRefreshSection()) {
    await refreshSection();
  }
});
els.mapToggle.addEventListener("change", () => {
  if (state.section.mode && els.mapToggle.checked) {
    els.mapToggle.checked = false;
    setStatus("Section picking uses Cartesian PPI view. Province overlay stayed off.", true);
  }
  refreshPlot();
});
els.prevButton.addEventListener("click", () => stepTimeline(-1));
els.nextButton.addEventListener("click", () => stepTimeline(1));
els.playButton.addEventListener("click", togglePlayback);
els.speedSelect.addEventListener("change", () => {
  if (state.playing) {
    togglePlayback();
    togglePlayback();
  }
});

els.sectionModeToggle.addEventListener("change", (event) => {
  if (event.target.checked) {
    beginSectionPicking(false);
    return;
  }
  setSectionMode(event.target.checked);
});
els.sectionPickButton.addEventListener("click", () => {
  beginSectionPicking(true);
});
els.sectionOpenButton.addEventListener("click", () => {
  if (state.section.objectUrl) {
    openSectionPanel(true);
  }
});
els.sectionClearInlineButton.addEventListener("click", () => {
  resetSectionSelection();
  setStatus("Section selection cleared.");
});
[
  els.sectionStartX,
  els.sectionStartY,
  els.sectionEndX,
  els.sectionEndY,
].forEach((input) => {
  input.addEventListener("change", readSectionInputs);
});
els.sectionRenderButton.addEventListener("click", () => {
  readSectionInputs();
  if (!sectionReady()) {
    setStatus("Set both section endpoints before rendering.", true);
    return;
  }
  refreshSection();
});
els.sectionPickAgainButton.addEventListener("click", () => {
  beginSectionPicking(true);
});
els.sectionClearButton.addEventListener("click", () => {
  resetSectionSelection();
  setStatus("Section selection cleared.");
});

els.canvas.addEventListener("wheel", (event) => {
  event.preventDefault();
  const factor = event.deltaY < 0 ? 1.12 : 0.9;
  applyZoom(event.clientX, event.clientY, factor);
}, { passive: false });

els.canvas.addEventListener("mousedown", (event) => {
  if (!state.canvas.image) {
    return;
  }
  state.canvas.dragging = true;
  state.canvas.dragMoved = false;
  state.canvas.dragStartX = event.clientX;
  state.canvas.dragStartY = event.clientY;
  state.canvas.dragOffsetX = state.canvas.offsetX;
  state.canvas.dragOffsetY = state.canvas.offsetY;
});

window.addEventListener("mousemove", (event) => {
  if (!state.canvas.dragging) {
    return;
  }
  if (Math.abs(event.clientX - state.canvas.dragStartX) > 3 || Math.abs(event.clientY - state.canvas.dragStartY) > 3) {
    state.canvas.dragMoved = true;
  }
  state.canvas.offsetX = state.canvas.dragOffsetX + (event.clientX - state.canvas.dragStartX);
  state.canvas.offsetY = state.canvas.dragOffsetY + (event.clientY - state.canvas.dragStartY);
  constrainCanvasView();
  drawCanvas();
});

window.addEventListener("mouseup", () => {
  state.canvas.dragging = false;
});

els.canvas.addEventListener("click", (event) => {
  if (!state.section.mode || !state.canvas.image || state.canvas.dragMoved) {
    return;
  }
  const point = clientToRadarKm(event.clientX, event.clientY);
  if (!point) {
    setStatus("Click inside the radar data square to place section endpoints.", true);
    return;
  }
  pickSectionPoint(point);
});

els.canvas.addEventListener("dblclick", () => {
  if (!state.canvas.image) {
    return;
  }
  resetCanvasView();
  drawCanvas();
});

window.addEventListener("resize", () => {
  if (state.canvas.image) {
    resetCanvasView();
    constrainCanvasView();
  }
  drawCanvas();
});

scanDirectory();

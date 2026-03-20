from easydict import EasyDict as edict

__C = edict()
cfg = __C

###set interp configure
__C.interp = edict()
__C.interp.mroi = 500 ##最小的影响半径
__C.interp.coeff = 50 ##随着距离圈增长，影响半径变大，增长倍数

###set 3D radar network configure
__C.network = edict()
__C.network.radar_dirs = []
__C.network.field_names = ["dBZ"]
__C.network.field_range_mode = {"dBZ": "native"}
__C.network.lon_min = None
__C.network.lon_max = None
__C.network.lat_min = None
__C.network.lat_max = None
__C.network.lon_res_deg = 0.01
__C.network.lat_res_deg = 0.01
__C.network.level_heights = []
__C.network.product_level_heights = []
__C.network.max_range_km = 460.0
__C.network.time_tolerance_minutes = 10
__C.network.composite_method = "quality_weighted"
__C.network.influence_radius_m = 300000.0
__C.network.quality_weight_terms = ["range", "beam", "vertical", "qc", "blockage"]
__C.network.range_decay_m = 230000.0
__C.network.beam_radius_ref_m = 3000.0
__C.network.vertical_gap_ref_deg = 0.5
__C.network.qc_weight_field = "QC_MASK"
__C.network.blockage_config = None
__C.network.fillvalue = -999.0
__C.network.output_dir = None
__C.network.effective_earth_radius = None
__C.network.file_pattern = "*.bin*"
__C.network.parallel = True
__C.network.max_workers = None
__C.network.blind_method = "hybrid"
__C.network.output_products = ["CR", "VIL", "ET"]
__C.network.plot_overview = True
__C.network.plot_style = "reference"
__C.network.plot_canvas_px = [920, 790]
__C.network.plot_height_levels = []
__C.network.plot_output_dir = None
__C.network.plot_format = "png"
__C.network.plot_product_for_levels = "CAPPI"
__C.network.coverage_policy = "observable_only"
__C.network.vil_min_dbz = 18.0
__C.network.vil_max_dbz_cap = 56.0
__C.network.et_threshold_dbz = 18.0
__C.network.use_qc = False
__C.network.qc_method = "dualpol"
__C.network.qc_band = "C"
__C.network.qc_use_existing_kdp = True
__C.network.qc_fallback = "original"
__C.network.qc_clear_air_mode = "label"
__C.network.qc_clear_air_max_ref = 15.0
__C.network.qc_clear_air_max_rhohv = 0.97
__C.network.qc_clear_air_max_phidp_texture = 10.0
__C.network.qc_clear_air_max_snr = 20.0

# -*- coding: utf-8 -*-
"""
Configuration file for the NuistRadar Object , modified from Py-art
注意：距离单位统一到m
"""

FILL_VALUE = -999.0
_LIGHT_SPEED = 2.99792458e8

DEFAULT_METADATA = {
    # Metadata for radar attributes. These closely follow the CF/Radial
    # standard

    'azimuth': {
        'units': 'degrees',
        'standard_name': 'beam_azimuth_angle',
        'long_name': 'azimuth_angle_from_true_north',
        'axis': 'radial_azimuth_coordinate',
        'comment': 'Azimuth of antenna relative to true north'},

    'elevation': {
        'units': 'degrees',
        'standard_name': 'beam_elevation_angle',
        'long_name': 'elevation_angle_from_horizontal_plane',
        'axis': 'radial_elevation_coordinate',
        'comment': 'Elevation of antenna relative to the horizontal plane'},

    'range': {
        'units': 'meters',
        'standard_name': 'projection_range_coordinate',
        'long_name': 'range_to_measurement_volume',
        'axis': 'radial_range_coordinate',
        'spacing_is_constant': 'true',
        'comment': (
            'Coordinate variable for range. Range to center of each bin.')},

    'time': {
        'standard_name': 'time',
        'long_name': 'time_in_seconds_since_volume_start',
        'calendar': 'gregorian',
        'comment': ('Coordinate variable for time. '
                    'Time at the center of each ray, in fractional seconds '
                    'since the global variable time_coverage_start')},

    'fixed_angle': {
        'long_name': 'Target angle for sweep',
        'units': 'degrees',
        'standard_name': 'target_fixed_angle'},

    'rays_per_sweep': {
        'long_name': 'Number of rays in each sweep',
        'units': 'count'},

    'latitude': {
        'long_name': 'Latitude',
        'standard_name': 'Latitude',
        'units': 'degrees_north'},

    'longitude': {
        'long_name': 'Longitude',
        'standard_name': 'Longitude',
        'units': 'degrees_east'},

    'altitude': {
        'long_name': 'Altitude',
        'standard_name': 'Altitude',
        'units': 'meters',
        'positive': 'up'},

    'nyquist_velocity': {
        'units': 'meters_per_second',
        'comments': "Unambiguous velocity",
        'meta_group': 'instrument_parameters',
        'long_name': 'Nyquist velocity'},

    'unambiguous_range': {
        'units': 'meters',
        'comments': 'Unambiguous range',
        'meta_group': 'instrument_parameters',
        'long_name': 'Unambiguous range'},
    'frequency': {
        'units': 'GHZ',
        'meta_group': 'instrument_parameters',
        'long_name': 'Radiation frequency'},

    # Reflectivity fields
    'reflectivity': {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'Reflectivity',
        'valid_max': 80.0,
        'valid_min': -30.0,
        'coordinates': 'elevation azimuth range'},
    "corrected_reflectivity": {
        'units': 'dBZ',
        'standard_name': 'corrected_equivalent_reflectivity_factor',
        'long_name': 'Corrected reflectivity',
        'valid_max': 80.0,
        'valid_min': -30.0,
        'coordinates': 'elevation azimuth range'},
    'scan_type':{
        'units': "string",
        'standard_name':"radar scan type",
        'long_name':"Type of scan, one of 'ppi', 'rhi', 'sector' or 'other'"
    },
    'total_power': {
        'units': 'dBZ',
        'standard_name': 'equivalent_reflectivity_factor',
        'long_name': 'Total power',
        'valid_max': 80.0,
        'valid_min': -30.0,
        'coordinates': 'elevation azimuth range'},

    # Velocity fields
    'velocity': {
        'units': 'meters_per_second',
        'standard_name': 'radial_velocity_of_scatterers_away_from_instrument',
        'long_name': 'Mean dopper velocity',
        'valid_max': 50.0,
        'valid_min': -50.0,
        'coordinates': 'elevation azimuth range'},

    "corrected_velocity": {
        'units': 'meters_per_second',
        'standard_name': 'corrected_radial_velocity_of_scatterers_away_from_instrument',
        'valid_max': 50.0,
        'valid_min': -50.0,
        'long_name': 'Corrected mean doppler velocity',
        'coordinates': 'elevation azimuth range'},

    # Spectrum width fields
    'spectrum_width': {
        'units': 'meters_per_second',
        'standard_name': 'doppler_spectrum_width',
        'long_name': 'Doppler spectrum width',
        'valid_max': 30.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    # Dual-polarization fields
    'differential_reflectivity': {
        'units': 'dB',
        'standard_name': 'log_differential_reflectivity_hv',
        'long_name': 'Differential reflectivity',
        'valid_max': 8.0,
        'valid_min': -2.0,
        'coordinates': 'elevation azimuth range'},

    'cross_correlation_ratio': {
        'units': 'ratio',
        'standard_name': 'cross_correlation_ratio_hv',
        'long_name': 'Cross correlation ratio (RHOHV)',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},

    'normalized_coherent_power': {
        'units': 'ratio',
        'standard_name': 'normalized_coherent_power',
        'long_name': 'Normalized coherent power',
        'valid_max': 1.0,
        'valid_min': 0.0,
        'comment': 'Also know as signal quality index (SQI)',
        'coordinates': 'elevation azimuth range'},

    'differential_phase': {
        'units': 'degrees',
        'standard_name': 'differential_phase_hv',
        'long_name': 'Differential phase (PhiDP)',
        'valid_max': 360.0,
        'valid_min': 0.0,
        'coordinates': 'elevation azimuth range'},
    'specific_differential_phase': {
        'units': 'degrees/km',
        'standard_name': 'specific_differential_phase_hv',
        'long_name': 'Specific differential phase (KDP)',
        'valid_max': 5,
        'valid_min': -2,
        'coordinates': 'elevation azimuth range'},
    'clutter_flag':{
        'standard_name': 'clutter_flag',
    },
    # Depolarization ratio fields
    'linear_depolarization_ratio': {
        'units': 'dB',
        'standard_name': 'log_linear_depolarization_ratio_hv',
        'long_name': 'Linear depolarization ratio',
        'valid_max': 0,
        'valid_min': -40.0,
        'coordinates': 'elevation azimuth range'},

}

# CINRAD files
CINRAD_field_mapping = {
    # moment: radar field name
    'dBT':'total_power',
    'dBZ': "reflectivity",
    'V': "velocity",
    'W': "spectrum_width",
    'SQI':'normalized_coherent_power',
    'CPA':None,
    'ZDR': "differential_reflectivity",
    'LDR': "linear_depolarization_ratio",
    'CC': "cross_correlation_ratio",
    'PhiDP': "differential_phase",
    'KDP': "specific_differential_phase",
    'CP':None,
    'FLAG':None,
    'HCL':None,
    'CF': "clutter_flag",
    'Zc': "corrected_reflectivity",
    'Vc': "corrected_velocity",
    'Wc': None,
}

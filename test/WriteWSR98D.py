import struct
import datetime
from netCDF4 import num2date

import numpy as np

from pycwr.io import read_auto
from pycwr.io.BaseDataProtocol.WSR98DProtocol import dtype_98D


# [dtype, scale, offset, binlength]
flag = {'total_power': [1, 2, 66, 1],
        'reflectivity': [2, 2, 66, 1],
        'velocity': [3, 2, 129, 1],
        'spectrum_width': [4, 2, 129, 1],
        'normalized_coherent_power': [5, 2, 129, 1],
        'CPA': [6, 2, 129, 1],
        'differential_reflectivity': [7, 16, 130, 1],
        'linear_polarization_ratio': [8, 16, 130, 1],
        'cross_correlation_ratio': [9, 200, 5, 1],
        'differential_phase': [10, 100, 50, 2],
        'specific_differential_phase': [11, 10, 50, 1],
        'CP': [12, 2, 129, 1],
        'radar_echo_classification': [14, 1, 0, 1],
        'CF': [15, 2, 129, 1],
        'signal_to_noise_ratio': [16, 2, 20, 1],
        'horizontal_signal_noise_ratio': [16, 2, 20, 1],
        'vertical_signal_noise_ratio': [17, 2, 20, 1],
        'corrected_reflectivity': [32, 2, 66, 1],
        'corrected_velocity': [33, 2, 129, 1],
        'corrected_spectrum_width': [34, 2, 129, 1],
        'corrected_differential_reflectivity': [35, 16, 130, 1]
        }

def date2julian(scan_start_time):
    """计算扫描开始时间
    return units:seconds
    """
    delta = scan_start_time - datetime.datetime(1970, 1, 1)
    return int(delta.total_seconds())

def date2julian_msec(scan_start_time):
    """计算扫描开始时间
    return units:seconds
    """
    delta = scan_start_time - datetime.datetime(1970, 1, 1)
    return int(delta.microseconds)

cut_type = np.dtype([('ProcessMode', '<i4'), ('WaveForm', '<i4'), ('PRF_1', '<f4'), ('PRF_2', '<f4'),
                     ('UnfoldMode', '<i4'), ('Azimuth', '<f4'), ('Elevation', '<f4'), ('StartAngle', '<f4'),
                     ('EndAngle', '<f4'), ('AngleResolution', '<f4'), ('ScanSpeed', '<f4'), ('LogResolution', '<i4'),
                     ('DopplerResolution', '<i4'), ('MaximumRange', '<i4'), ('MaximumRange2', '<i4'), ('StartRange', '<i4'),
                     ('Sample_1', '<i4'), ('Sample_2', '<i4'), ('PhaseMode', '<i4'), ('AtmosphericLoss', '<f4'),
                     ('NyquistSpeed', '<f4'), ('MomentsMask', '<i8'), ('MomentsSizeMask', '<i8'), ('MiscFilterMask', '<i4'),
                     ('SQIThreshold', '<f4'), ('SIGThreshold', '<f4'), ('CSRThreshold', '<f4'), ('LOGThreshold', '<f4'),
                     ('CPAThreshold', '<f4'), ('PMIThreshold', '<f4'), ('DPLOGThreshold', '<f4'), ('ThresholdsReserved', 'V4'),
                     ('dBTMask', '<i4'), ('dBZMask', '<i4'), ('Velocity', '<i4'), ('SpectrumWidthMask', '<i4'), ('ZDRMask', '<i4'),
                     ('MaskResvered', 'V12'), ('ScanSync', '<i4'), ('Direction', '<i4'), ('GroundClutterClassifierType', '<i2'),
                     ('GroundClutterFilterType', '<i2'), ('GroundClutterFilterNotchWidth', '<i2'),
                     ('GroundClutterFilterWindow', '<i2'), ('Reserved', 'V72')])

def write_CN98D(ArtRd, filename, SiteID = b"Z0000", SiteName = b"CMA", BandwithH = 0.93, BandwithV = 0.93, RDAVersion = 722951,
                RadarType = 1, TaskName = b"VCP21", freq = 2.75, nyquist_vel = 27.8):

    sweep_end_ray_index = ArtRd.sweep_end_ray_index["data"]
    sweep_start_ray_index = ArtRd.sweep_start_ray_index["data"]
    rays_per_sweep = sweep_end_ray_index - sweep_start_ray_index + 1
    vals = list(ArtRd.fields.keys())

    for i, ival in enumerate(vals):
        if not ival in flag.keys():
            vals.pop(ival)

    elevation = ArtRd.elevation["data"]
    azimuth = ArtRd.azimuth["data"]
    ranges = ArtRd.range["data"]
    fields = ArtRd.fields
    fix_angle = ArtRd.fixed_angle["data"]
    if "frequency" in ArtRd.instrument_parameters.keys():
        Frequency = ArtRd.instrument_parameters["frequency"]["data"][0]/10**9
    else:
        Frequency = freq

    if "nyquist_velocity" in ArtRd.instrument_parameters.keys():
        nyquist_velocity = ArtRd.instrument_parameters["nyquist_velocity"]["data"]
    else:
        nyquist_velocity = np.array([nyquist_vel,] * len(rays_per_sweep),)

    latitude = ArtRd.latitude["data"][0]
    longitude = ArtRd.longitude["data"][0]
    altitude = ArtRd.altitude["data"][0]
    scan_time = np.asarray(num2date(ArtRd.time["data"], ArtRd.time["units"]))
    StartTime = date2julian(scan_time[0])
    res_bins = int(ranges[1] -ranges[0])

    fmt_GenericHeader = '<' + ''.join([i[1] for i in dtype_98D.BaseDataHeader['GenericHeaderBlock']])
    fmt_SiteConfiguration = '<' + ''.join([i[1] for i in dtype_98D.BaseDataHeader['SiteConfigurationBlock']])
    fmt_TaskConfiguration = '<' + ''.join([i[1] for i in dtype_98D.BaseDataHeader['TaskConfigurationBlock']])
    fmt_RadialHeader = '<' + ''.join([i[1] for i in dtype_98D.RadialHeader()])
    fmt_RadialData = '<' + ''.join([i[1] for i in dtype_98D.RadialData()])

    buffer=b''
    GenericHeader = struct.pack(fmt_GenericHeader, 1297371986, -1, -1, 1, -1, b"")
    buffer += GenericHeader
    SiteConfiguration = struct.pack(fmt_SiteConfiguration, b"%s" % SiteID, b"%s" % SiteName, latitude, longitude,
                                    int(altitude), int(altitude), Frequency, BandwithH, BandwithV, RDAVersion,
                                    RadarType, b"")
    buffer += SiteConfiguration
    TaskConfiguration = struct.pack(fmt_TaskConfiguration, b"%s" % TaskName, b"", -1, -1, -1, StartTime,
                                    len(rays_per_sweep), -1, -1, -1, -1, -1, -1, -1, -1, -1, b"")
    buffer += TaskConfiguration

    CutConfig = []
    for isweep in range(len(rays_per_sweep)):
        CutConfig.append(np.array([(-1, -1, -1., -1., -1, -1., fix_angle[isweep], -1., -1., BandwithH, -1., res_bins, res_bins,
                                    int(ranges[-1]), -2147483648, 0, 28, -2147483648, 0, 0.011, nyquist_velocity[isweep], 69254, 1024, 510, 0.4, 5.,
                                    60., 3., 0., 0.45, 5., b'', 1, 1, 1, 1, 32, b'', 0, 1, 3, 1, 3, 1, b'')], dtype=cut_type))

    CutConfiguration = np.concatenate(CutConfig).tobytes()
    buffer += CutConfiguration

    count = 0
    for isweep in range(len(rays_per_sweep)):
        for irays in range(rays_per_sweep[isweep]):
            if count == 0:
                RadialState = 3
            elif (irays == (rays_per_sweep[isweep] - 1)) and (isweep == (len(rays_per_sweep) - 1)):
                RadialState = 4
            elif irays == (rays_per_sweep[isweep] - 1):
                RadialState = 2
            elif irays == 0:
                RadialState = 0
            else:
                RadialState = 1
            sum_bytes = 0
            for ivar in vals:
                sum_bytes += flag[ivar][3] * len(ranges) + 32

            tmp_radial_header = struct.pack(fmt_RadialHeader, RadialState, 0, count + 1, irays + 1, isweep + 1,
                                            azimuth[count], elevation[count], int(date2julian(scan_time[count])),
                                            date2julian_msec(scan_time[count]), sum_bytes, len(vals), b"")
            buffer += tmp_radial_header

            for ivar in vals:
                tmp_radial_data = struct.pack(fmt_RadialData, flag[ivar][0], flag[ivar][1], flag[ivar][2],
                                                flag[ivar][3], -32768, flag[ivar][3] * len(ranges), b"")
                buffer += tmp_radial_data
                encode = (np.asarray(fields[ivar]["data"][count]) * flag[ivar][1]) + flag[ivar][2]

                if flag[ivar][3] == 1:
                    code = np.where(np.isnan(encode), 3, encode).astype("u1")
                elif flag[ivar][3] == 2:
                    code = np.where(np.isnan(encode), 3, encode).astype("u2")
                else:
                    raise ValueError("Bin Length Must Be 1 or 2!")
                tmp_data = code.tobytes()
                buffer += tmp_data
            count += 1
    with open(filename, "wb") as f:
        f.write(buffer)
    return

if __name__ == "__main__":
    import glob
    import os, time
    import pyart

    files = glob.glob("/Volumes/LYUFC-DRIVE/NEXRAD_Files/*_V06")
    for fileread in files:
        s = time.time()
        filesave = os.path.join("/Volumes/LYUFC-DRIVE/outfiles/", os.path.basename(fileread) + ".bin")
        # fileread = r"/Users/zhengyu/Downloads/Z_RADR_I_Z9515_20160623063100_O_DOR_SA_CAP.bin"
        # basedata = WSR98DFile.WSR98DBaseData(filename)
        ArtRd = pyart.io.read_nexrad_archive(fileread)
        # filesave = r"/Users/zhengyu/Downloads/test.bin"
        write_CN98D(ArtRd, filesave)
        e = time.time()
        print(f"{fileread} runing {(e - s)} seconds" )
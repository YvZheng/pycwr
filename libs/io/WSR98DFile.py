# -*- coding: utf-8 -*-
import sys

sys.path.append("../../")
import struct
import numpy as np
from BaseDataProtocol.WSR98DProtocol import dtype_98D
from util import _prepare_for_read, _unpack_from_buf, julian2date_SEC
import time
from libs.core.NRadar import NuistRadar

class WSR98DBaseData(object):
    """
    解码新一代双偏振的数据格式
    """

    def __init__(self, filename):
        super(WSR98DBaseData, self).__init__()
        self.fid = _prepare_for_read(filename)  ##对压缩的文件进行解码
        self._check_standard_basedata()  ##确定文件是standard文件
        self.header = self._parse_BaseDataHeader()
        self.radial = self._parse_radial()
        self.nrays = len(self.radial)
        status = np.array([istatus['RadialState'] for istatus in self.radial[:]])
        self.sweep_start_ray_index = np.where((status == 0) | (status == 3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.nsweeps = len(self.sweep_start_ray_index)
        self.fid.close()

    def _check_standard_basedata(self):
        """
        :param fid: file fid
        :return:
        """
        assert self.fid.read(4) == b'RSTM', 'file in not a stardand WSR-98D file!'
        self.fid.seek(0, 0)
        return

    def _parse_BaseDataHeader(self):
        BaseDataHeader = {}
        fixed_buf = self.fid.read(dtype_98D.CutConfigurationBlockPos)  ##读取前面固定头的信息

        BaseDataHeader['GenericHeader'], _ = _unpack_from_buf(fixed_buf, \
                                                              dtype_98D.GenericHeaderBlockPos,
                                                              dtype_98D.BaseDataHeader['GenericHeaderBlock'])
        BaseDataHeader['SiteConfig'], _ = _unpack_from_buf(fixed_buf, \
                                                           dtype_98D.SiteConfigurationBlockPos,
                                                           dtype_98D.BaseDataHeader['SiteConfigurationBlock'])
        BaseDataHeader['TaskConfig'], _ = _unpack_from_buf(fixed_buf,
                                                           dtype_98D.TaskConfigurationBlockPos,
                                                           dtype_98D.BaseDataHeader['TaskConfigurationBlock'])
        cut_buf = self.fid.read(dtype_98D.CutConfigurationBlockSize * \
                                BaseDataHeader['TaskConfig']['CutNumber'])
        BaseDataHeader['CutConfig'] = np.frombuffer(cut_buf, dtype_98D.BaseDataHeader['CutConfigurationBlock'])
        return BaseDataHeader

    def _parse_radial(self):
        radial = []
        buf = self.fid.read(dtype_98D.RadialHeaderBlockSize)
        while len(buf) == dtype_98D.RadialHeaderBlockSize:  ##read until EOF
            RadialDict, _ = _unpack_from_buf(buf, 0, dtype_98D.RadialHeader())
            self.MomentNum = RadialDict['MomentNumber']
            self.LengthOfData = RadialDict['LengthOfData']
            RadialDict['fields'] = self._parse_radial_single()
            radial.append(RadialDict)
            buf = self.fid.read(dtype_98D.RadialHeaderBlockSize)
        return radial

    def _parse_radial_single(self):
        radial_var = {}
        for _ in range(self.MomentNum):
            Mom_buf = self.fid.read(dtype_98D.MomentHeaderBlockSize)
            Momheader, _ = _unpack_from_buf(Mom_buf, 0, dtype_98D.RadialData())
            Data_buf = self.fid.read(Momheader['Length'])
            assert (Momheader['BinLength'] == 1) | (Momheader['BinLength'] == 2), "Bin Length has problem!"
            if Momheader['BinLength'] == 1:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u1", offset=dtype_98D.MomentHeaderBlockSize)).astype(int)
            else:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u2", offset=dtype_98D.MomentHeaderBlockSize)).astype(int)
            radial_var[dtype_98D.flag2Product[Momheader['DataType']]] = np.where(dat_tmp >= 5, \
                                                                                 (dat_tmp - Momheader['Offset']) /
                                                                                 Momheader['Scale'],
                                                                                 np.nan).astype(np.float32)
        return radial_var

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return np.concatenate([[nyquist, ] * ray for nyquist, ray in zip(self.header['CutConfig']['NyquistSpeed'], \
                                                                         self.get_rays_per_sweep())], axis=0)

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate([[nyquist, ] * ray for nyquist, ray in zip(self.header['CutConfig']['MaximumRange'], \
                                                                         self.get_rays_per_sweep())], axis=0)

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return np.array([julian2date_SEC(iray['Seconds'], iray['MicroSeconds']) for iray in self.radial])

    def get_sweep_end_ray_index(self):
        """
        获取每个sweep的结束的index，包含在内
        :return:(nsweep)
        """
        return self.sweep_end_ray_index

    def get_sweep_start_ray_index(self):
        """
        获取每个sweep的开始的index
        :return:(nsweep)
        """
        return self.sweep_start_ray_index

    def get_rays_per_sweep(self):
        """
        获取每个sweep的径向数
        :return:(nsweep)
        """
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.array([self.radial[iray]['Azimuth'] for iray in range(self.nrays)])

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return np.array([self.radial[iray]['Elevation'] for iray in range(self.nrays)])

    def get_latitude_longitude_altitude_frequency(self):
        """
        获取经纬度高度，雷达频率
        :return:lat, lon, alt, frequency
        """
        return self.header['SiteConfig']['Latitude'], self.header['SiteConfig']['Longitude'], \
               self.header['SiteConfig']['Height'], self.header['SiteConfig']['Frequency'] / 1000.

    def get_scan_type(self):
        if self.header['TaskConfig']['ScanType'] in [0, 1]:
            return "ppi"
        elif self.header['TaskConfig']['ScanType'] in [2, 5]:
            return "rhi"
        elif self.header['TaskConfig']['ScanType'] in [3, 4]:
            return "sector"
        else:
            return "other"


class WSR98D2NRadar(object):
    """到NusitRadar object 的桥梁"""

    def __init__(self, WSR98D):
        self.WSR98D = WSR98D
        assert np.all(self.WSR98D.header['CutConfig']['LogResolution'] == \
                      self.WSR98D.header['CutConfig']['DopplerResolution']), "dop not match dBZ!"
        v_index_alone = self.get_v_idx()
        dBZ_index_alone = self.get_dbz_idx()
        assert np.all(v_index_alone == (dBZ_index_alone + 1)), """v and dBZ not equal!"""
        for index_with_dbz, index_with_v in zip(dBZ_index_alone, v_index_alone):
            assert (self.WSR98D.header["CutConfig"]["Elevation"][index_with_v] - \
                    self.WSR98D.header["CutConfig"]["Elevation"][index_with_dbz]) < 0.5, \
                "warning! maybe it is a problem."
            self.interp_dBZ(index_with_dbz, index_with_v)
        ind_remove = self.get_reomve_radial_num()
        self.radial = [iray for ind, iray in enumerate(self.WSR98D.radial) if ind not in ind_remove]
        self.nsweeps = self.WSR98D.nsweeps - dBZ_index_alone.size
        status = np.array([istatus['RadialState'] for istatus in self.radial[:]])
        self.sweep_start_ray_index = np.where((status == 0) | (status == 3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.nrays = len(self.radial)
        self.scan_type = self.WSR98D.get_scan_type()
        self.latitude, self.longitude, self.altitude, self.frequency = \
            self.WSR98D.get_latitude_longitude_altitude_frequency()
        self.header = self.WSR98D.header
        self.header['CutConfig'] = np.delete(self.header['CutConfig'], dBZ_index_alone)
        self.bins_per_sweep = self.get_nbins_per_sweep()
        self.max_bins = self.bins_per_sweep.max()
        self.range = self.get_range_per_radial(self.max_bins)
        self.azimuth = self.get_azimuth()
        self.elevation = self.get_elevation()
        self.fields = self._get_fields()

    def get_reomve_radial_num(self):
        """获得需要remove的radial的index"""
        dBZ_alone = self.get_dbz_idx()
        index_romove = []
        for isweep in dBZ_alone:
            index_romove.extend(range(self.WSR98D.sweep_start_ray_index[isweep], \
                                      self.WSR98D.sweep_end_ray_index[isweep] + 1))
        return index_romove

    def get_v_idx(self):
        """获取需要插值的sweep, 插值到有径向速度仰角"""
        flag = np.array([(("V" in self.WSR98D.radial[idx]['fields'].keys()) and (
                    "dBZ" not in self.WSR98D.radial[idx]['fields'].keys())) \
                         for idx in self.WSR98D.sweep_start_ray_index])
        return np.where(flag == 1)[0]

    def get_dbz_idx(self):
        """获取含有dbz的sweep"""
        flag = np.array([(("dBZ" in self.WSR98D.radial[idx]['fields'].keys()) and (
                    "V" not in self.WSR98D.radial[idx]['fields'].keys())) \
                         for idx in self.WSR98D.sweep_start_ray_index])
        return np.where(flag == 1)[0]

    def interp_dBZ(self, field_with_dBZ_num, field_without_dBZ_num):
        """
        将dBZ插值到不含dBZ的仰角
        :param field_with_dBZ_num: 要插值的sweep num, （从0开始）
        :param field_without_dBZ_num: 要插值到的sweep num, (从0开始)  which to evaluate the interpolated values
        :return:
        """
        azimuth = self.WSR98D.get_azimuth() ##获取98D的方位角
        assert (field_with_dBZ_num + 1) == field_without_dBZ_num, "check interp sweep!"
        dbz_az = azimuth[self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]: \
                         self.WSR98D.sweep_end_ray_index[field_with_dBZ_num] + 1]
        v_az = azimuth[self.WSR98D.sweep_start_ray_index[field_without_dBZ_num]: \
                       self.WSR98D.sweep_end_ray_index[field_without_dBZ_num] + 1]
        dbz_idx = np.argmin(np.abs(dbz_az.reshape(-1, 1) - v_az.reshape(1, -1)), axis=0) + \
                  self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]

        v_idx = np.arange(self.WSR98D.sweep_start_ray_index[field_without_dBZ_num], \
                          self.WSR98D.sweep_end_ray_index[field_without_dBZ_num] + 1)
        keys = self.WSR98D.radial[self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]]['fields'].keys()
        for ind_dbz, ind_v in zip(dbz_idx, v_idx):
            for ikey in keys:
                self.WSR98D.radial[ind_v]["fields"][ikey] = self.WSR98D.radial[ind_dbz]["fields"][ikey]

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.array([self.radial[iray]['Azimuth'] for iray in range(self.nrays)])

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return np.array([self.radial[iray]['Elevation'] for iray in range(self.nrays)])

    def get_rays_per_sweep(self):
        """
        获取每个sweep的径向数
        :return:(nsweep)
        """
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return np.array([julian2date_SEC(iray['Seconds'], iray['MicroSeconds']) for iray in self.radial])

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return np.concatenate([[nyquist, ] * ray for nyquist, ray in zip(self.header['CutConfig']['NyquistSpeed'], \
                                                                         self.get_rays_per_sweep())], axis=0)

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate([[nyquist, ] * ray for nyquist, ray in zip(self.header['CutConfig']['MaximumRange'], \
                                                                         self.get_rays_per_sweep())], axis=0)

    def get_sweep_end_ray_index(self):
        """
        获取每个sweep的结束的index，包含在内
        :return:(nsweep)
        """
        return self.sweep_end_ray_index

    def get_sweep_start_ray_index(self):
        """
        获取每个sweep的开始的index
        :return:(nsweep)
        """
        return self.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """
        确定每个sweep V探测的库数
        :return:
        """
        return np.array([self.radial[idx]['fields']['V'].size for idx in self.sweep_start_ray_index])

    def get_range_per_radial(self, length):
        """
        确定径向每个库的距离 range变量
        :param length:
        :return:
        """
        Resolution = self.header['CutConfig']['DopplerResolution'][0]
        return np.linspace(Resolution, Resolution * length, length)

    def _get_fields(self):
        """将所有的field的数据提取出来"""
        fields = {}
        field_keys = self.radial[0]['fields'].keys()
        for ikey in field_keys:
            fields[ikey] = np.array([self._add_or_del_field(iray['fields'], ikey) for iray in self.radial])
        return fields

    def _add_or_del_field(self, dat_fields, key):
        """
        根据fields的key提取数据
        :param dat_fields: fields的数据
        :param key: key words
        :return:
        """
        length = self.max_bins
        if key not in dat_fields.keys():
            return np.full((length,), np.nan)
        dat_ray = dat_fields[key]
        assert dat_ray.ndim == 1, "check dat_ray"
        if dat_ray.size >= length:
            return dat_ray[:length]
        else:
            out = np.full((length,), np.nan)
            out[:dat_ray.size] = dat_ray
            return out

    def get_NRadar_nyquist_speed(self):
        """array shape (nsweeps)"""
        return self.header['CutConfig']['NyquistSpeed']

    def get_NRadar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.header['CutConfig']['MaximumRange']

    def get_fixed_angle(self):
        if self.scan_type == "rhi":
            return self.header['CutConfig']['Azimuth']
        else:
            return self.header['CutConfig']['Elevation']

    def ToNuistRadar(self):
        """将WSR98D数据转为Nuist Radar的数据格式"""
        return NuistRadar(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_NRadar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_NRadar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps)

    def ToPyartRadar(self):
        pass


if __name__ == "__main__":
    start = time.time()
    # test = WSR98DBaseData(r"E:\RadarBaseData\郑玉\利奇马\LQM_Z9513_NT\Z_RADR_I_Z9513_20190811234540_O_DOR_SAD_CAP_FMT.bin.bz2")
    test = WSR98DBaseData(r"E:\RadarBaseData\StandardFormat\厦门\Z9592.20160728.111443.AR2.bz2")
    # test = WSR98DBaseData(r'E:\RadarBaseData\StandardFormat\NUIST\NUIST.20170408.153319.AR2')
    wsr98d = WSR98D2NRadar(test)
    NRadar = wsr98d.ToNuistRadar()
    end = time.time()
    print(end - start)

# -*- coding: utf-8 -*-
import struct
import numpy as np
from BaseDataProtocol.StandardProtocol import dtype_standard
from util import _prepare_for_read, _unpack_from_buf, julian2date_SEC
import time


class WSR98DBaseData(object):
    """
    解码新一代双偏振的数据格式
    """

    def __init__(self, filename):
        super(StandardBaseData, self).__init__()
        self.fid = _prepare_for_read(filename)  ##对压缩的文件进行解码
        self._check_standard_basedata()  ##确定文件是standard文件
        self.header = self._parse_BaseDataHeader()
        self.radial = self._parse_radial()
        self.nrays = len(self.radial)
        status = np.array([istatus['RadialState'] for istatus in self.radial[:]])
        self.sweep_start_ray_index = np.where((status == 0) | (status == 3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.nseeps = len(self.sweep_start_ray_index)
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
        fixed_buf = self.fid.read(dtype_standard.CutConfigurationBlockPos)  ##读取前面固定头的信息

        BaseDataHeader['GenericHeader'], _ = _unpack_from_buf(fixed_buf, \
                                                              dtype_standard.GenericHeaderBlockPos,
                                                              dtype_standard.BaseDataHeader['GenericHeaderBlock'])
        BaseDataHeader['SiteConfig'], _ = _unpack_from_buf(fixed_buf, \
                                                           dtype_standard.SiteConfigurationBlockPos,
                                                           dtype_standard.BaseDataHeader['SiteConfigurationBlock'])
        BaseDataHeader['TaskConfig'], _ = _unpack_from_buf(fixed_buf,
                                                           dtype_standard.TaskConfigurationBlockPos,
                                                           dtype_standard.BaseDataHeader['TaskConfigurationBlock'])
        cut_buf = self.fid.read(dtype_standard.CutConfigurationBlockSize * \
                                BaseDataHeader['TaskConfig']['CutNumber'])
        BaseDataHeader['CutConfig'] = np.frombuffer(cut_buf, dtype_standard.BaseDataHeader['CutConfigurationBlock'])
        return BaseDataHeader

    def _parse_radial(self):
        radial = []
        buf = self.fid.read(dtype_standard.RadialHeaderBlockSize)
        while len(buf) == dtype_standard.RadialHeaderBlockSize:  ##read until EOF
            RadialDict, _ = _unpack_from_buf(buf, 0, dtype_standard.RadialHeader())
            self.MomentNum = RadialDict['MomentNumber']
            self.LengthOfData = RadialDict['LengthOfData']
            RadialDict['fields'] = self._parse_radial_single()
            radial.append(RadialDict)
            buf = self.fid.read(dtype_standard.RadialHeaderBlockSize)
        return radial

    def _parse_radial_single(self):
        radial_var = {}
        for _ in range(self.MomentNum):
            Mom_buf = self.fid.read(dtype_standard.MomentHeaderBlockSize)
            Momheader, _ = _unpack_from_buf(Mom_buf, 0, dtype_standard.RadialData())
            Data_buf = self.fid.read(Momheader['Length'])
            assert (Momheader['BinLength'] == 1) | (Momheader['BinLength'] == 2), "Bin Length has problem!"
            if Momheader['BinLength'] == 1:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u1", offset=dtype_standard.MomentHeaderBlockSize)).astype(int)
            else:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u2", offset=dtype_standard.MomentHeaderBlockSize)).astype(int)
            radial_var[dtype_standard.flag2Product[Momheader['DataType']]] = np.where(dat_tmp >= 5, \
                                                                                      (dat_tmp - Momheader['Offset']) /
                                                                                      Momheader['Scale'],
                                                                                      np.nan).astype(np.float32)
        return radial_var

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return np.concatenate([[nyquist,] * ray for nyquist, ray in zip(self.header['CutConfig']['NyquistSpeed'],\
                                                                 self.get_rays_per_sweep())], axis=0)

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate([[nyquist,] * ray for nyquist, ray in zip(self.header['CutConfig']['MaximumRange'],\
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
               self.header['SiteConfig']['Height'], self.header['SiteConfig']['Frequency']/1000.


class WSR98D2NRadar(object):
    """到NusitRadar object 的桥梁"""
    def __init__(self, filename):

        pass


if __name__ == "__main__":
    start = time.time()
    test = StandardBaseData(r"E:\RadarBaseData\郑玉\利奇马\LQM_Z9513_NT\Z_RADR_I_Z9513_20190809012256_O_DOR_SAD_CAP_FMT.bin.bz2")
    #test = StandardBaseData(r"E:\RadarBaseData\StandardFormat\厦门\Z9592.20160728.111443.AR2.bz2")
    end = time.time()
    print(end - start)

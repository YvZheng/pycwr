# -*- coding: utf-8 -*-
import struct
import numpy as np
from BaseDataProtocol.SABProtocol import dtype_sab
from util import _prepare_for_read, _unpack_from_buf, julian2date
import time

class SABBaseData(object):
    """
    解码SA/SB/CB/SC2.0的雷达数据，仅仅对数据（dBZ, V, W）做了转换
    """

    def __init__(self, filename):
        super(SABBaseData, self).__init__()
        self.fid = _prepare_for_read(filename)
        self.RadialNum, self.nrays = self._RadialNum_SAB_CB() ##检查文件有无问题
        self.radial = self._parse_radial()
        status = np.array([istatus['RadialStatus'] for istatus in self.radial[:]])
        self.sweep_start_ray_index = np.where((status==0)|(status==3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.nsweeps = len(self.sweep_start_ray_index)
        self.fid.close()

    def _RadialNum_SAB_CB(self):
        """f: a file-like object was provided, 确定雷达数据的径向字节长度"""
        assert self.fid.read(28)[14:16] == b'\x01\x00', 'file in not a valid SA/SB/CB file!'
        self.fid.seek(0, 0) ##移动到开头
        data_len = len(self.fid.read())
        assert (data_len%2432 == 0) | (data_len%4132 == 0) |(data_len%3132 == 0), "file size has problems!"
        ###判断雷达数据类型SA/SB 或者 CB
        if data_len % 2432 == 0:
            RadialNum = 2432
            self.Type = "SAB"
        elif data_len%4132 == 0:
            RadialNum = 4132
            self.Type = 'CB'
        else:
            RadialNum = 3132
            self.Type = 'SC'
        self.fid.seek(0, 0) ##移动到开头
        return RadialNum, int(data_len/RadialNum)

    def _parse_radial(self):
        """
        循环读取所有径向数据
        :param fid:
        :return:
        """
        radial = []
        print(self.nrays)
        for _ in range(self.nrays):
            radial.append(self._parse_radial_single(self.fid.read(self.RadialNum)))
        return radial

    def _parse_radial_single(self, radial_buf):
        Radial = {}
        RadialHeader, size_tmp = _unpack_from_buf(radial_buf, 0, dtype_sab.RadialHeader())
        Radial.update(RadialHeader)
        RadialDataDtype = dtype_sab.RadialData(RadialHeader['GatesNumberOfReflectivity'],
                                               RadialHeader['GatesNumberOfDoppler'])
        FieldSize = RadialHeader['GatesNumberOfReflectivity'] + RadialHeader['GatesNumberOfDoppler']*2
        RadialData = np.frombuffer(radial_buf[size_tmp:size_tmp+FieldSize], dtype=RadialDataDtype)
        Radial['fields'] = {}
        Radial['fields']['dBZ'] = np.where(RadialData['dBZ']>1, (RadialData['dBZ'].astype(int) - 2)/2.-32,
                                           np.nan).astype(np.float32)
        Radial['fields']['V'] = np.where(RadialData['V'] > 1, (RadialData['V'].astype(int) - 2) / 2. - 63.5,
                                           np.nan).astype(np.float32)
        Radial['fields']['W'] = np.where(RadialData['W'] > 1, (RadialData['W'].astype(int) - 2) / 2. - 63.5,
                                         np.nan).astype(np.float32)
        return Radial

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return np.array([iradial['Nyquist']/100. for iradial in self.radial])
    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离 units:km
        :return:(nRays)
        """
        return np.array([iradial['URange']/10. for iradial in self.radial])
    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return np.array([julian2date(iradial['JulianDate'], iradial['mSends']) for iradial in self.radial])
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

    def get_bins_per_sweep(self):
        """
        获取每个sweep的库数
        :return: REF(nsweep), DOP(nsweep)
        """
        REF = np.array([(self.radial[iray]['fields']['dBZ']).size for iray in self.sweep_start_ray_index])
        DOP = np.array([(self.radial[iray]['fields']['V']).size for iray in self.sweep_start_ray_index])
        return REF, DOP

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.array([iradial['AZ']/8.* 180./4096. for iradial in self.radial])
    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return np.array([iradial['El']/8.* 180./4096. for iradial in self.radial])

    def get_latitude_longitude_altitude_frequency(self):
        """
        获取经纬度高度，雷达频率
        :return:lat, lon, alt, frequency
        """
        pass

class Standard2NRadar(object):
    """到NusitRadar object 的桥梁"""
    def __init__(self):
        pass



if __name__ == "__main__":
    start = time.time()
    #test = SABBaseData("/home/zy/data/code_data/ERIC/Radar/2010081202.13A")
    test = SABBaseData("/home/zy/Desktop/Z9240/Z_RADR_I_Z9240_20190703080123_O_DOR_SC_CAP.bin")
    end = time.time()
    print(end-start)
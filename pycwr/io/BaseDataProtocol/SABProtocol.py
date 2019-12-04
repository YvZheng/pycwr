import numpy as np

UINT8 = "B"
UINT16 = "H"
UINT32 = "I"

class SABFormat(object):

    def __init__(self):
        super(SABFormat, self).__init__()
        self.RadialHeaderSize = 128
        self.InfSize = 28

    def RadialHeader(self):
        return (
            ('reserve0', '14s'),
            ('flag', UINT16), ##1-标识雷达数据
            ('reserve1', '12s'),
            ('mSends', UINT32),  # 径向数据收集时间(毫秒,自 00:00 开始)
            ('JulianDate', UINT16),  # 儒略日(Julian)表示,自 1970 年 1 月 1 日开始
            ('URange', UINT16),  # 不模糊距离(表示:数值/10.=千米)
            ('AZ', UINT16),  # 方位角(编码方式:[数值/8.]*[180./4096.]=度)
            ('RadialNumber', UINT16),  # 当前仰角内径向数据序号
            ('RadialStatus', UINT16),  # 径向数据状态
            ('El', UINT16),
            # 仰角 (编码方式:[数值/8.]*[180./4096.]=度)
            ('ElNumber', UINT16),  # 体扫内的仰角数
            ('RangeToFirstGateOfRef', UINT16),  # 反射率数据的第一个距离库的实际距离(单位:米)
            ('RangeToFirstGateOfDop', UINT16),  # 多普勒数据的第一个距离库的实际距离(单位:米)
            ('GateSizeOfReflectivity', UINT16),  # 反射率数据的距离库长(单位:米)
            ('GateSizeOfDoppler', UINT16),  # 多普勒数据的距离库长(单位:米)
            ('GatesNumberOfReflectivity', UINT16),  # 反射率的距离库数
            ('GatesNumberOfDoppler', UINT16),  # 多普勒的距离库数
            ('CutSectorNumber', UINT16),  # 扇区号
            ('CalibrationConst', UINT32),  # 系统订正常数
            ('PtrOfReflectivity', UINT16),  # 反射率数据指针(偏离雷达数据信息头的字节数) 表示第一个反射率数据的位置
            ('PtrOfVelocity', UINT16),  # 速度数据指针(偏离雷达数据信息头的字节数),表示第一个速度数据的位置
            ('PtrOfSpectrumWidth', UINT16),  # 谱宽数据指针(偏离雷达数据信息头的字节数),表示第一个谱宽数据的位置
            ('ResolutionOfVelocity', UINT16),  # 多普勒速度分辨率。 2:表示 0.5 米/秒
            ('VcpNumber', UINT16),  # 体扫(VCP)模式
            ('reserve2', '14s'),  # 保留
            ('Nyquist', UINT16),  # Nyquist 速度(表示:数值/100. = 米/秒)
            ('reserve3', '38s'))

    def RadialData(self, rnumber, vnumber):
        """
        :param rnumber: 反射率的库数
        :param vnumber: 多普勒的库数
        :param radialbins: 整个径向的长度
        :return:
        """
        return np.dtype([('dBZ', 'u1', rnumber),
                ('V', 'u1', vnumber),
                ('W', 'u1', vnumber)])

dtype_sab = SABFormat()
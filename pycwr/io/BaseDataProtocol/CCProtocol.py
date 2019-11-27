# -*- coding: utf-8 -*-
import numpy as np

INT8 = 'b'
SHORT = 'h'
INT = 'i'
FLOAT = 'f'

UINT8 = "B"
UINT16 = "H"
UINT32 = "I"


class CCFormat(object):

    def __init__(self):
        super(CCFormat, self).__init__()
        self.PerRadialSize = 3000  # 每根径向的字节长度
        self.BaseDataHeaderSize = 1024
        self.HeaderSize1 = 218
        self.HeaderSize1_pos = 0
        self.CutSize = 660
        self.CutSize_pos = 218
        self.HeaderSize2 = 146
        self.HeaderSize2_pos = 878

        self.BaseDataHeader = dict(
            RadarHeader1=(
            ('cFileType', '16s'),             ##3830数据标识(CINRADC)
            ('cCountry', '30s'),              ##国家名
            ('cProvince', '20s'),             ##省名
            ('cStation', '40s'),              ##站名
            ('cStationNumber', '10s'),        ##区站号
            ('cRadarType', '20s'),            ##雷达型号
            ('cLongitude', '16s'),            ##天线所在经度
            ('cLatitude', '16s'),             ##天线所在纬度
            ('lLongitudeValue', INT),         ##具体经度
            ('lLatitudeValue', INT),          ##具体纬度
            ('lHeight', INT),                 ##天线海拔高度
            ('sMaxAngle', SHORT),             ##地物阻挡最大仰角
            ('sOptAngle', SHORT),             ##最佳观测仰角
            ('ucSYear1', UINT8),              ##观测开始时间的年千百位(19-20)
            ('ucSYear2', UINT8),              ##观测开始时间的年十个位(00-99)
            ('ucSMonth', UINT8),              ##观测开始时间的月(1-12)
            ('ucSDay', UINT8),                ##观测开始时间的日(1-31)
            ('ucSHour', UINT8),               ##观测开始时间的时(0-23)
            ('ucSMinute', UINT8),             ##观测开始时间的分(0-59)
            ('ucSSecond', UINT8),             ##观测开始时间的秒(0-59)
            ('ucTimeFrom', UINT8),            ##时间来源 0-计算机时钟(1天内未对时)
                                              ##		   1-计算机时钟(1天内已对时)
            ('ucEYear1', UINT8),              ##		   2-GPS
            ('ucEYear2', UINT8),              ##		   3-其它
            ('ucEMonth', UINT8),              ##观测结束时间的年千百位(19-20)
            ('ucEDay', UINT8),                ##观测结束时间的年十个位(00-99)
            ('ucEHour', UINT8),               ##观测结束时间的月(1-12)
            ('ucEMinute', UINT8),             ##观测结束时间的日(1-31)
            ('ucESecond', UINT8),             ##观测结束时间的时(0-23)
            ('ucScanMode', UINT8),            ##观测结束时间的分(0-59)
                                              ##观测结束时间的秒(0-59)
            ('ulSmilliSecond', UINT32),       ##扫描方式  1-RHI
            ('usRHIA', UINT16),               ##		   10-PPI和ZPPI
                                              ##		   1XX=VPPI(XX为扫描圈数)
            ('sRHIL', SHORT),                 ##以微秒为单位表示的秒的小数位
                                              ##RHI所在的方位角(0.01度为单位)
            ('sRHIH', SHORT),                 ## PPI和VPPI时为FFFF
                                              ##RHI所在的最低仰角(0.01度为单位)
            ('usEchoType', UINT16),           ##PPI和VPPI时为FFFF
                                              ##RHI所在的最高仰角(0.01度为单位)
            ('usProdCode', UINT16),           ##PPI和VPPI时为FFFF
                                              ##回波类型  0x405a-Z  0x406a-V  0x407a-W
            ('ucCalibration', UINT8),         ##		   0x408a-ZVW三要素
                                              ##数据类型  0x8001-PPI数据  0x8002-RHI数据
            ('remain1', '3s'),                ##   0x8003-VPPI数据  0x8004-单强度PPI数据
        ),
        CutConfigX30 = np.dtype([
            ('usMaxV', "u2"),                 ##最大可测速度(厘米/秒)
            ('usMaxL', "u2"),                 ##最大可测距离(10米)
            ('usBindWidth', "u2"),            ##库长(米)
            ('usBinNumber', "u2"),            ##每径向库数
            ('usRecordNumber', "u2"),         ##本圈径向数
            ('usArotate', "u2"),              ##本圈转速(0.01度/秒)
            ('usPrf1', "u2"),                 ##本圈第一次重复频率(0.1Hz)对应单重频或双重频的高者
            ('usPrf2', "u2"),                 ##本圈第二次重复频率(0.1Hz)对应双重频的低者
            ('usSpulseW', "u2"),              ##本圈脉宽(微秒)
            ('usAngle', "i2"),                ##仰角(0.01度)
            ('cSweepStatus', "u1"),           ##1=单要素	2=三要素(单重频)	3=三要素(双重频)
            ('cAmbiguousp', "u1"),            ##0=无软件退模糊	1=软件退模糊
        ]),
        RadarHeader2 = (
            ("remain2", "2s"),                ##保留字节
            ('lAntennaG', INT),               ##		   3-1月内人工
            ('lPower', INT),                  ##保留字
            ('lWavelength', INT),             ##保留字, 放VPPISCANPARAMETER数据
            ('usBeamH', UINT16),              ##该结构的说明见后
            ('usBeamL', UINT16),              ##天线增益(0.001dB)
            ('usPolarization', UINT16),       ##峰值功率(瓦)
                                              ##波长(微米)
            ('usLogA', UINT16),               ##垂直波束宽度(秒)
            ('usLineA', UINT16),              ##水平波束宽度(秒)
            ('usAGCP', UINT16),               ##极化状态 0-水平 1-垂直 2-双偏振
            ('usFreqMode', UINT16),           ## 		  3-圆偏振 4-其它
                                              ##对数动态范围(0.01dB)
            ('usFreqRepeat', UINT16),         ##线性动态范围(0.01dB)
            ('usPPPPulse', UINT16),           ##AGC延迟量(微秒)
            ('usFFTPoint', UINT16),           ##频率方式	1-单重复频率  2-双重复频率3:2
            ('usProcessType', UINT16),        ##			3-双重复频率4:3
                                              ##重复频率
            ('ucClutterT', UINT8),            ##PPP脉冲数
            ('cSidelobe', INT8),              ##FFT间隔点数
            ('ucVelocityT', UINT8),           ##信号处理方式	1-PPP	2-全程FFT
            ('ucFilderP', UINT8),             ##				3-单库FFT
                                              ##杂波消除阀值(即STC)
            ('ucNoiseT', UINT8),              ##第一旁瓣(dB)
            ('ucSQIT', UINT8),                ##速度门限
            ('ucIntensityC', UINT8),          ##地物消除方式	0-无		1-IIR滤波器1
                                              ##			2-IIR滤波器2	3-IIR滤波器3
            ('ucIntensityR', UINT8),          ##			4-IIR滤波器4
                                              ##噪声消除阀值(即强度门限)
            ('ucCalNoise', UINT8),            ##SQI门限
            ('ucCalPower', UINT8),            ##DVIP强度值估算采用的通道
            ('ucCalPulseWidth', UINT8),       ##	 1-对数通道 2-线性通道
            ('ucCalWorkFreq', UINT8),         ##强度值估算是否距离订正
            ('ucCalLog', UINT8),              ## 0-无(dB) 1-已订正(dBZ)
            ('remain3', '92s'),               ##噪声系数标定值
            ('liDataOffset', UINT32),         ##发射功率标定值
            ('remain4', "1s"),                ##保留字
        ))
    def RadialData(self, radialnumber):
        """
        :param radialnumber: 库数
        :return:
        """
        return np.dtype([('dBZ', 'i2', radialnumber),
                ('V', 'i2', radialnumber),
                ('W', 'i2', radialnumber)])

dtype_cc = CCFormat()



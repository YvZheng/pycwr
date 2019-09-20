import numpy as np

INT8 = "b"
INT16 = "h"
INT32 = "i"
UINT8 = "B"
UINT16 = "H"
UINT32 = "I"

class SCFormat(object):

    def __init__(self):
        super(SCFormat, self).__init__()
        self.BaseDataHeaderSize = 1024##头的字节长度
        self.PerRadialSize = 4000#每根径向的字节长度
        self.RadarSitePos = 0
        self.RadarPerformanceParamPos = 170
        self.RadarObserationParamPos_1 = 201
        self.LayerParamPos = 217
        self.RadarObserationParamPos_2 = 847
        self.BaseDataHeader = dict(
        RadarSite = (
        ('country', '30s'),   #国家名，文本格式输入
        ('province', '20s'),  #省名，文本格式输入
        ('station', '40s'),   #站名，文本格式输入
        ('stationnumber', '10s'), #区站号，文本格式输入
        ('radartype', '20s'), #雷达型号，文本格式输入
        ('longitude', '16s'), #天线所在经度，文本格式输入.书写格式例：E 115°32′12″
        ('latitude', '16s'),  #天线所在纬度，文本格式输入.书写格式例：N 35°30′15″
        ('longitudevalue', INT32), #天线所在经度的数值，以毫秒为计数单位（十进制）东经（E）为正，西经（W）为负 ！！后续转换
        ('latitudevalue', INT32),  #天线所在纬度的数值，以毫秒为计数单位（十进制）北纬（N）为正，南纬（S）为负 ！！后续转换
        ('height', INT32),         #天线的海拔高度以毫米为计数单位   ！！后续转换
        ('Maxangle', INT16),       #测站四周地物阻挡的最大仰角（以秒为计数单位） ！！后续转换
        ('Opangle', INT16),        #测站的最佳观测仰角（地物回波强度<10dbz，以秒为计数单位） ！！后续转换
        ('MangFreq', INT16)       #磁控管频率（通过此频率可计算雷达波长）//Whyan
        ),
        RadarPerformanceParam = (
        ('AntennaG', INT32),         #天线增益，以0.001db为计数单位
        ('BeamH', UINT16),           #垂直波束宽度，以微秒为计数单位
        ('BeamL', UINT16),           #水平波束宽度，以微秒为计数单位
        ('polarizations',UINT8),     #极化状况 # 0 = 水平 1 = 垂直 2 = 双偏振 3 = 圆偏振 4 = 其他
        ('sidelobe', INT8),          #第一旁瓣计数单位：db（注意：输入负号）
        ('Power', UINT32),           #雷达脉冲峰值功率，以瓦为计数单位
        ('wavelength',UINT32),       #波长，以微米为计数单位
        ('logA', UINT16),            #对数接收机动态范围,以0.01db为计数单位
        ('LineA', UINT16),           #线性接收机动态范围,以0.01为计数单位
        ('AGCP', UINT16),            #AGC延迟量，以微秒为计数单位
        ('clutterT',UINT8),          #杂波消除阀值，计数单位0.01db
        ('VelocityP',UINT8),         #速度处理方式 0 = 无速度处理 1 = PPP 2 = FFT 3 = 随机编码
        ('filderP',UINT8),           #地物消除方式 0 = 无地物消除 1 = 地物杂波图扣除法 2 = 地物杂波图 + 滤波器处理 3 = 滤波器处理 4 = 谱分析处理
        ('noiseT',UINT8),            #噪声消除阀值	（0-255）
        ('SQIT',UINT8),              #SQI阀值，以0.01为计数单位
        ('intensityC',UINT8),        #rvp强度值估算采用的通道 1 = 对数通道 2 = 线性通道
        ('intensityR',UINT8),        #强度估算是否进行了距离订正 0 = 无 1 = 已进行了距离订正
        ),
        RadarObserationParam_1 = (
        ('stype',UINT8),#  扫描方式
                        #  1 = RHI
                        #  10 = PPI
                        #1XX = Vol 	XX为扫描圈数
        ('syear', UINT16),  #观测记录开始时间的年的十位个位（01-99）
        ('smonth' , UINT8), #观测记录开始时间的月（1-12）
        ('sday'   , UINT8), #观测记录开始时间的日（1-31）
        ('shour'  , UINT8), #观测记录开始时间的时（00-23）
        ('sminute', UINT8), #观测记录开始时间的分（00-59）
        ('ssecond', UINT8), #观测记录开始时间的秒（00-59）
        ('Timep'  , UINT8), #时间来源
                            #0 = 计算机时钟，但一天内未进行对时
                            #1 = 计算机时钟，但一天内已进行对时
                            #2 = GPS
                            #3 = 其他
        ('smillisecond', UINT32), #秒的小数位（计数单位微秒）
        ('calibration', UINT8),   #	标校状态
                                  #	 0 = 无标校
                                  #	 1 = 自动标校
                                  #	 2 = 1星期内人工标校
                                  #	 3 = 1月内人工标校
                                  #	其他码不用
        ('intensityI', UINT8), #强度积分次数（32-128）
        ('VelocityP', UINT8), ),#速度处理样本数（31-255）(样本数-1）
        LayerParamX30 = np.dtype([
        ('ambiguousp', 'u1'), #本层退模糊状态  !!这里重复30次
                               # 0 = 无退模糊状态
                               # 1 = 软件退模糊
                               # 2 = 双T退模糊
                               # 3 = 批式退模糊
                               # 4 = 双T + 软件退模糊
                               # 5 = 批式 + 软件退模糊
                               # 6 = 双PPI退模糊
                               # 9 = 其他方式
        ('Arotate', 'u2'),    #本层天线转速,计数单位:0.01度/秒
        ('Prf1', 'u2'),       #本层的第一种脉冲重复频率,计数单位: 1/10 Hz
        ('Prf2', 'u2'),       #	本层的第二种脉冲重复频率,计数单位: 1/10 Hz
        ('spulseW', 'u2'),    #	本层的脉冲宽度,计数单位:	微秒
        ('MaxV', 'u2'),       #	本层的最大可测速度,计数单位:	厘米/秒
        ('MaxL', 'u2'),        #本层的最大可测距离，以10米为计数单位
        ('binWidth', 'u2'),    #本层数据的库长，以分米为计数单位
        ('binnumber', 'u2'),   #本层扫描线水平方向的点数
        ('recordnumber', 'u2'),#  本层扫描线垂直方向的点数 	本层径向数
        ('Swangles', 'i2'),     #  本层的仰角，计数单位	：1/100度
        ]),
        RadarObserationParam_2 = (
        ('RHIA', UINT16),  #作RHI时的所在方位角，计数单位为1/100度, 作其他扫描时不用
        ('RHIL', INT16),   #作RHI时的最低仰角，计数单位为1/100度 , 作其他扫描时不用
        ('RHIH', INT16),   #作RHI时的最高仰角，计数单位为1/100度, 作其他扫描时不用
        ('Eyear',UINT16), #	观测结束时间的年的十位个位（01-99）
        ('Emonth',UINT8), #	观测结束时间的月（1-12）
        ('Eday',UINT8),   #	观测结束时间的日（1-31）
        ('Ehour',UINT8),  #	观测结束时间的时（00-23）
        ('Eminute',UINT8),#	观测结束时间的分（00-59）
        ('Esecond',UINT8),#	观测结束时间的秒（00-59）
        ('Etenth',UINT8), ),)#	观测结束时间的1/100秒（00-59） cFileType[16];

    def RadialHeader(self):

        return (('sStrAz', UINT16), #方位角计算方法： 实际Az= Az*360./65536
        ('sStrEl', UINT16), #仰角计算方法： 实际El= El*120./65536
        ('sEndAz', UINT16),
        ('sEndEl', UINT16),)

    def RadialData(self):
        """每个径向的数据提取"""
        return np.dtype(
        [('dBZ', 'u1'),
        ('V', 'u1'),
        ('dBT','u1'),
        ('W', 'u1')])

dtype_sc = SCFormat()

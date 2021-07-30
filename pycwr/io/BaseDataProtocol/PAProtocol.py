# -*- coding: utf-8 -*-
import numpy as np

SHORT = 'h'
INT = 'i'
FLOAT = 'f'
LONG = "q"
UINT16 = "H"
UINT32 = "I"

class PAFormat(object):

    def __init__(self):
        self.GenericHeaderBlockPos = 0
        self.SiteConfigurationBlockPos = 32
        self.TaskConfigurationBlockPos = 160
        self.BeamConfigurationBlockPos = 416
        self.BeamConfigurationBlockSize = 384
        self.CutConfigurationBlockSize = 256
        self.RadialHeaderBlockSize = 128
        self.MomentHeaderBlockSize = 32
        self.flag2Product = {1:'dBT', 2:'dBZ', 3:'V', 4:'W', 5:'SQI', 6:'CPA', 7:'ZDR', 8:'LDR',\
                             9:'CC', 10:'PhiDP', 11:'KDP', 12:'CP', 13:'FLAG', 14:'HCL', 15:'CF',\
                             16:'SNRH', 17:'SNRV', 18:'Flag', 19:"Flag", 20:'Flag',21:"Flag",\
                             22:'Flag', 23:"Flag", 24:'Flag', 25:'Flag',26:"Flag", 27:'Flag', \
                             28:'Flag', 29:"Flag", 30:'Flag',31:'Flag', 32:"Zc", 33:'Vc', 34:'Wc',\
                             35:'ZDRc', 0:'Flag'}
        self.BaseDataHeader = dict(
        GenericHeaderBlock = (
            ('MagicWord', INT),  # 固定标志，用来指示雷达数据文件。
            ('MajorVersion', SHORT),  # 主版本号
            ('MinorVersion', SHORT),  # 次版本号
            ('GenericType', INT),  # 文件类型 1-基数据文件 2-气象产品文件
            ('ProductType', INT),  # 产品类型 文件类型1时无效
            ('Reserved01', '16s')  # 保留字段
        ),
        SiteConfigurationBlock = (
            ('SiteCode', '8s'),  # 站号具有唯一性，用来区别不同的雷达
            ('SiteName', '32s'),  # 站点名称，如NUIST
            ('Latitude', FLOAT),  # 雷达站天线所在位置纬度 degree
            ('Longitude', FLOAT),  # 雷达站天线所在位置经度 degree
            ('Height', INT),  # 天线馈源水平时海拔高度  meter
            ('Ground', INT),  # 雷达塔楼地面海拔高度  meter
            ('Frequency', FLOAT),  # 工作频率  Mhz
            ('BeamWidthHori', FLOAT),  # 水平波束宽度 degree
            ('BeamWidthVert', FLOAT),  # 垂直波束宽度 degree
            ('RDAVersion', INT),  # RDA版本号
            ('RadarType', SHORT),  # 雷达类型 SA/SB/SC/SAD
            ('Reserved02', '54s')  # 保留字段
        ),
        TaskConfigurationBlock = (
            ('TaskName', '32s'),  # 任务名称，如VCP21
            ('TaskDescription', '128s'),  # 任务描述
            ('PolarizationType', INT),  # 极化方式 1-水平 2-垂直 3-水平/垂直同时 4-水平/垂直交替
            ('ScanType', INT),  # 扫描任务类型 0 – 体扫 1–单层PPI 2 – 单层RHI 3 – 单层扇扫 4 – 扇体扫 5 – 多层RHI 6 – 手工扫描
            ('BeamNumber', INT),  # M
            ('CutNumber', INT),  # N
            ('RayOrder', INT),  # 扫描层数 根据扫描任务类型确定的扫描层数
            ('VolumeStartTime', LONG),  # 扫描开始时间 扫描开始时间为UTC标准时间计数,1970年1月1日0时为起始计数基准点 units:seconds
            ('Reserved03', '68s'),  # 保留字段
        ),
        BeamConfigurationBlock = np.dtype([
            ('BeamIndex', 'i4'),  # 基数据仰角编号
            ('BeamType', 'i4'),  # 对应的发射波束索引
            ('SubPulseNumber', 'i4'),  # PPI模式的俯仰角
            ('TxBeamDirection', 'f4'),  # 在本接收方向上，发射波束的增益
            ('TxBeamWidthH', 'f4'),  #
            ('TxBeamWidthV', 'f4'),  #
            ('TxBeamGain', 'f4'),  #
            ('Reserved00', '100V'),  # 处理模式     1-PPP 2-FFT
            ('SubPulseStrategy', 'i4'),  # 基数据仰角编号
            ('SubPulseModulation', 'i4'),  # 对应的发射波束索引
            ('SubPulseFrequency', 'f4'),  # PPI模式的俯仰角
            ('SubPulseBandWidth', 'f4'),  # 在本接收方向上，发射波束的增益
            ('SubPulseWidth', 'i4'),  #
            ('Reserved01', '236V'),  # 处理模式     1-PPP 2-FFT
            ]),
        CutConfigurationBlock = np.dtype([
            ('CutIndex', 'i2'),  # 基数据仰角编号
            ('TxBeamIndex', 'i2'),  # 对应的发射波束索引
            ('Elevation', 'f4'),  # PPI模式的俯仰角
            ('TxBeamGain', 'f4'),  # 在本接收方向上，发射波束的增益
            ('RxBeamWidthH', 'f4'),  #
            ('RxBeamWidthV', 'f4'),  #
            ('RxBeamGain', 'f4'),  #
            ('ProcessMode', 'i4'),  # 处理模式     1-PPP 2-FFT
            ('WaveForm', 'i4'),  # 波形类别
            ('N1_PRF_1', 'f4'),  # 脉冲重复频率1 双PRF表示高PRF 单PRF表示唯一值 Hz
            ('N1_PRF_2', 'f4'),  # 脉冲重复频率2 双PRF表示低PRF 单PRF无效 Hz
            ('N2_PRF_1', 'f4'),  # 脉冲重复频率1 双PRF表示高PRF 单PRF表示唯一值 Hz
            ('N2_PRF_2', 'f4'),  # 脉冲重复频率2 双PRF表示低PRF 单PRF无效 Hz
            ('UnfoldMode', 'i4'),  # 速度退模糊方法1 – 单PRF 2 –双PRF3:2模式 3 –双PRF4:3模式 4 –双PRF 5:4模式
            ('Azimuth', 'f4'),  # 方位角 degree
            ('StartAngle', 'f4'),  # 起始角度 degree
            ('EndAngle', 'f4'),  # 起始角度 degree
            ('AngleResolution', 'f4'),  # 角度分辨率 degree
            ('ScanSpeed', 'f4'),  # 扫描速度   degree/sec
            ('LogResolution', 'f4'),  # 强度分辨率 强度数据的距离分辨率 meter
            ('DopplerResolution', 'f4'),  # 多普勒数据的距离分辨率 meter
            ('MaximumRange', 'i4'),  # 对应PRF1的最大探测距离 meter
            ('MaximumRange2', 'i4'),  # 对应PRF2的最大探测距离  meter
            ('StartRange', 'i4'),  # 数据探测的起始距离 meter
            ('Sample_1', 'i4'),  # 对应于脉冲重复频率1的采样个数
            ('Sample_2', 'i4'),  # 对应于脉冲重复频率2的采样个数
            ('PhaseMode', 'i4'),  # 相位编码模式1 – 固定相位2 – 随机相位 3 – SZ编码
            ('AtmosphericLoss', 'f4'),  # 双程大气衰减值，精度为小数点后保留6位 dB/km
            ('NyquistSpeed', 'f4'),  # 理论最大不模糊速度  m/s
            ('MomentsMask', 'i8'),  # 数据类型掩码 0–不允许获取数据1 –允许获取数据。
            ('MomentsSizeMask', 'i8'),  # 数据大小掩码
            ('MiscFilterMask', 'i4'),  # 滤波设置掩码
            ('SQIThreshold', 'f4'),  # SQI门限
            ('SIGThreshold', 'f4'),  # SIG门限
            ('CSRThreshold', 'f4'),  # CSR门限
            ('LOGThreshold', 'f4'),  # LOG门限
            ('CPAThreshold', 'f4'),  # CPA门限
            ('PMIThreshold', 'f4'),  # PMI门限
            ('DPLOGThreshold', 'f4'),  # DPLOG门限
            ('ThresholdsReserved', '4V'),  # 保留字段
            ('dBTMask', 'i4'),  # dBT质控掩码 0应用 1未应用
            ('dBZMask', 'i4'),  # dBZ质控掩码
            ('Velocity', 'i4'),  # VEL质控掩码
            ('SpectrumWidthMask', 'i4'),  # SW质控掩码
            ('ZDRMask', 'i4'),  # ZDR质控掩码
            ('MaskResvered', '12V'),  # 保留
            ('ScanSync', '4V'),  # 扫描同步标志   用于多部雷达同步
            ('Direction', 'i4'),  # 天线运行方向 ppi有效 0顺时针 1逆时针
            ('GroundClutterClassifierType', 'i2'),  # 地物杂波图类型
            ('GroundClutterFilterType', 'i2'),  # 地物滤波类型
            ('GroundClutterFilterNotchWidth', 'i2'),  # 地物滤波宽度
            ('GroundClutterFilterWindow', 'i2'),  # 滤波窗口类型
            ('Reserved', '44V')  # 保留
        ])
        )
    def RadialHeader(self):
        RadialHeaderBlock = (
            ('RadialState', INT),  # 径向数据状态
            ('SpotBlank', INT),  # 消隐标志  0正常 1消隐
            ('SequenceNumber', INT),  # 序号，每个体扫径向从1开始计数
            ('RadialNumber', INT),  # 径向数，每个扫描从1开始
            ('ElevationNumber', INT),  # 仰角编号，从1开始
            ('Azimuth', FLOAT),  # 方位角  degree
            ('Elevation', FLOAT),  # 仰角   degree
            ('Seconds', LONG),  # 径向数据采集时间 UTC 从19700101 00：00开始 秒
            ('MicroSeconds', INT),  # 径向数据采集时间除去秒后留下的毫秒数
            ('LengthOfData', INT),  # 本径向数据块的长度
            ('MomentNumber', INT),  # 径向数据类别数量
            ('ScanBeamIndex', SHORT),  # 波束编号
            ('HorizontalEstimatedNoise', SHORT),  # 径向的水平通道估计噪声
            ('VerticalEstimatedNoise', SHORT),  # 径向数据类别数量
            ('PRFFLAG', UINT32),  # 波束编号
            ('Reserved04', '70s'),  # 保留
        )
        return RadialHeaderBlock

    def RadialData(self):
        MomentHeaderBlock = (
            ('DataType', INT),  # 数据类型
            ('Scale', INT),  # 比例
            ('Offset', INT),  # 偏移
            ('BinLength', SHORT),  # 单个库字节长度
            ('Flags', SHORT),  # 数据标志位，暂不使用
            ('Length', INT),  # 距离库的长度，不包括头 bytes
            ('Reserved05', '12s')  # 保留
        )
        return MomentHeaderBlock

dtype_PA = PAFormat()
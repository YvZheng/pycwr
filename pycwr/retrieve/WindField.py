import numpy as np

def VVP(azimuth, elevation, PPI_Vr, az_num, bin_num, fillvalue=-999.):
    """
    VVP方法反演风场，基于均匀风假设
    :param azimuth:np.ndarray,1d, 一维的方位角，对应PPI_Vr的第一维度
    :param elevation:constant, 常量，该层PPI扫描的仰角
    :param PPI_Vr:PPI扫描的径向速度，已退模糊
    :param az_num:反演取样的体积nazs
    :param bin_num:反演取样的体积nbins
    :param fillvalue:径向速度的缺测值
    :return:
    """
    assert az_num % 2 == 1, "az_num应为奇数"
    assert bin_num % 2 == 1, "bin_num应为奇数"

    Naz, NRange = PPI_Vr.shape
    half_bins = int(bin_num / 2)
    half_azs = int(az_num / 2)
    CosAz = np.cos(np.deg2rad(azimuth)).reshape(-1, 1)
    SinAz = np.sin(np.deg2rad(azimuth)).reshape(-1, 1)
    CosE = np.cos(np.deg2rad(elevation))

    D11 = np.pad((SinAz * CosE * SinAz * CosE).repeat(NRange, axis=1), ((half_azs, half_azs), (0, 0)), mode="wrap")
    D12 = np.pad((CosAz * CosE * SinAz * CosE).repeat(NRange, axis=1), ((half_azs, half_azs), (0, 0)), mode="wrap")
    D21 = np.pad((SinAz * CosE * CosAz * CosE).repeat(NRange, axis=1), ((half_azs, half_azs), (0, 0)), mode="wrap")
    D22 = np.pad((CosAz * CosE * CosAz * CosE).repeat(NRange, axis=1), ((half_azs, half_azs), (0, 0)), mode="wrap")

    C1 = np.pad(PPI_Vr * SinAz * CosE, ((half_azs, half_azs), (0, 0)), mode="wrap")
    C2 = np.pad(PPI_Vr * CosAz * CosE, ((half_azs, half_azs), (0, 0)), mode="wrap")
    Vr = np.pad(PPI_Vr, ((half_azs, half_azs), (0, 0)), mode="wrap")

    wind_u = np.full_like(PPI_Vr, fillvalue, dtype=np.float32)
    wind_v = np.full_like(PPI_Vr, fillvalue, dtype=np.float32)

    A = np.zeros((2, 2))
    B = np.zeros(2)
    for i in range(half_azs, Naz + half_azs):
        for j in range(half_bins, NRange - half_bins):
            flag = Vr[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1] != fillvalue
            if np.sum(flag) > 0.8 * az_num * bin_num:
                d11 = D11[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1]
                d12 = D12[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1]
                d21 = D21[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1]
                d22 = D22[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1]
                c1 = C1[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1]
                c2 = C2[i - half_azs:i + half_azs + 1, j - half_bins:j + half_bins + 1]
                A[0, 0] = np.sum(d11[flag])
                A[0, 1] = np.sum(d12[flag])
                A[1, 0] = np.sum(d21[flag])
                A[1, 1] = np.sum(d22[flag])
                B[0] = np.sum(c1[flag])
                B[1] = np.sum(c2[flag])
                x = np.linalg.solve(A, B)
                wind_u[i - half_azs, j] = x[0]
                wind_v[i - half_azs, j] = x[1]
    return wind_u, wind_v
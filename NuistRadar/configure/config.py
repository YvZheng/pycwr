from easydict import EasyDict as edict

__C = edict()
cfg = __C

###set interp configure
__C.interp = edict()
__C.interp.mroi = 500 ##最小的影响半径
__C.interp.coeff = 50 ##随着距离圈增长，影响半径变大，增长倍数

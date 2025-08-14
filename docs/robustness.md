# 鲁棒性与错误学（WSR98D 第一阶段）

本页描述 I/O 解析层的防御式策略、异常体系与可配置开关（默认保持旧行为）。

## 异常类型
- `CINRADFormatError`: 头部魔数/结构非法、关键字段不一致（如 `BinLength` 非 1/2，或 V/dBZ 扫描不匹配且未启用修复）。
- 后续将补充：`MissingSweepError`, `MissingRayError`, `EndianMismatchError`, `RecordLengthError`, `CalibrationMismatchError`。

## 开关参数
- `strict`: bool，默认 False。对潜在问题（如仰角差异≥0.5°）从告警提升为异常。
- `repair_missing`: bool，默认 False。对缺失/不一致的索引与字段进行可选修复：
  - 对 V/dBZ 扫描不对齐：插值 dBZ 至 V 方位角；无法配对的 sweep 以 NaN 填充（后续在 Py-ART 字段成为 Mask）。
  - 对 start/end 索引不等长：截取共同长度并保证 `end >= start`。

API 示例：
```python
from pycwr.io import read_WSR98D, read_auto
prd = read_WSR98D(path, strict=False, repair_missing=True)
prd2 = read_auto(path, strict=True, repair_missing=False)
```

## 防御式解析
- 头部魔数严格校验：非 `b'RSTM'` → `CINRADFormatError`。
- Moment `BinLength` 仅允许 1/2：否则 `CINRADFormatError`。
- V/dBZ 对齐：非 VCP26D/27D 任务，若 V 与 dBZ 非相邻 sweep：
  - `repair_missing=False` → `CINRADFormatError`
  - `repair_missing=True` → 插值/填充 NaN

## 掩码策略与数值一致性
- 修复分支对缺口填充 NaN；导出 Py-ART 字段为 `numpy.ma.MaskedArray`（NaN 转为 mask），保证下游统计稳定。
- 默认参数保持旧版行为不变；开启 `strict/repair_missing` 时的差异需在回归测试中给出黄金断言与误差阈值。


# 鲁棒性与错误学（占位）

本页将在 io-robustness PR 中补充：
- 防御式解析策略（魔数/字节序/记录长度）。
- 异常类型：`CINRADFormatError`, `MissingSweepError`, `MissingRayError`, `EndianMismatchError`, `RecordLengthError`, `CalibrationMismatchError`。
- 开关参数：`robust`, `fill_missing`, `strict_header` 与默认行为。
- 掩码策略与质量标志。


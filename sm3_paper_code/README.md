# SM3 论文门级代码

此目录是 `main-jcr-demo.tex` 当前资源统计对应的独立代码副本。它只包含实现和验证所需的 Python 文件，不包含早期统计报告、临时输出或绘图记录。原 `resource_estimation` 目录保留供追溯。

## 范围与统计口径

- 实现一个 512 比特消息分组的 64 轮压缩线路，以及两条状态线路的交替计算；包含 52 次 `W16` 至 `W67` 消息扩展、每轮 `W'` 的生成与还原、P0/P1、QFA、FF/GG、异地加法器的正算与反算，以及四个并行原地加法器。
- 直接统计发出的门；对图 3(b) 型 QAND† 的条件修正采用固定预算。加法器内部 QAND† 的条件修正沿用论文表 2、表 3 的口径，不另加到门数中。
- 输入消息和初始状态按符号寄存器分配；不计入其初值制备所需的 X 门。每轮 T 常量的设置和清理由程序实际计入。
- 不包含消息填充与分组结束后的最终前态异或。论文的 2176 层 T 深度是调度估算，程序目前没有独立计算该深度。

在 ProjectQ 0.8.0 环境中运行完整资源统计：

```powershell
python .\sm3_full_gate_level.py --rounds 64
```

当前输出应为：

| 指标 | 数量 |
| --- | ---: |
| T/T† 门 | 629248 |
| CNOT 门 | 1682448 |
| 单比特 Clifford 门 | 674304 |
| 测量 | 157312 |
| 峰值量子比特 | 2024 |

这些数值对应论文表 1 和第 4 节当前的门数及量子比特总数。论文的量子比特**分项时序说明**仍与代码的实际峰值时序不同；代码峰值为 `512` 位消息、`512` 位两条状态线路、`64` 位 T 寄存器、`256` 位当前 COMP 输出、`64` 位临时复制结果，以及 `616` 位四个原地加法器工作空间。

## 文件

| 文件 | 用途 |
| --- | --- |
| `sm3_full_gate_level.py` | 64 轮门级连接与资源统计入口 |
| `optimal_depth_adder.py`、`parallel_in_place_adders.py` | 原地加法器和四路并行工作空间 |
| `reversible_out_of_place.py` | 异地加法器正算与反算 |
| `sm3_gate_blocks.py` | QFA、FF、GG 及其反算 |
| `linear_permutations.py`、`message_schedule.py` | 68-CNOT P0/P1 网络与 16 字循环消息扩展 |
| `fixed_branch_accounting.py` | QAND† 固定条件修正预算 |
| `verify_*.py` | 门级部件、异地加法器与整体布尔功能的验证 |

可分别执行验证：

```powershell
python .\linear_permutations.py
python .\message_schedule.py
python .\verify_gate_blocks.py
python .\verify_reversible_adder.py
python .\verify_full_classical.py --rounds 3
```

`verify_full_classical.py` 将 QAND 临时替换为逻辑 Toffoli，只检查布尔接线与反算；资源数以 `sm3_full_gate_level.py` 为准。

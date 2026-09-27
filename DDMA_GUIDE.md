# DDMA 教学仿真：虚构配置与原始 Bin 顺序

入口为 [simulate_ddma.py](simulate_ddma.py)，使用 [radar_cfg_example.h](configs/radar_cfg_example.h) 中完全虚构的波形和相位码。模型是 4TX 同时发射、1RX 等效接收。

```powershell
py -3.11 simulate_ddma.py
```

配图、检测坐标和复数通道输出见 [ddma_results/REPORT.md](ddma_results/REPORT.md)。

## 1. 相位编码和频移

DDMA 为各 TX 设置逐 chirp 的线性相位变化，使回波产生不同的慢时间频移。[MathWorks DDMA 示例](https://www.mathworks.com/help/radar/ug/simulate-an-automotive-4d-imaging-mimo-radar.html) 介绍了这一调制方式和未占用频带的作用。

样例假设 256 个相位码对应一圈。若初相位码为 a、步进为 b：

\[
\phi_q[m]=2\pi(a_q+m b_q)/256.
\]

采用 `tx_ref * conj(rx)` 混频，编码在基带中为 `exp(-j phi_q[m])`。对于默认 **256 点** Doppler FFT：

| TX | 初相位码 | 初相位 | 步进码 | 步进角度 | 编码搬移（mod 256） |
|---|---|---:|---|---:|---:|
| TX1 | 0x10 | 22.5° | 0x00 | 0° | 0 |
| TX2 | 0x50 | 112.5° | 0x30 | 67.5° | 208 |
| TX3 | 0xA0 | 225° | 0x70 | 157.5° | 144 |
| TX4 | 0xE0 | 315° | 0xD0 | 292.5° | 48 |

这些码是独立构造的教学值，周期为 16 个 chirp。初相位影响复数通道相位，步进决定频谱搬移。

## 2. 二维 FFT 不重新排序

`process_ddma()` 保留 FFT 原始输出顺序，不调用 `fftshift`。二维图横轴为 **Range Bin**，纵轴为 **Doppler Bin**。

默认显示范围：Range Bin 0～127，Doppler Bin 0～255。Doppler Bin 0 为 DC，1～127 为正频率，128～255 按 FFT 约定对应负频率。

解码后的 k 用 `k_signed = k if k < N/2 else k-N` 换算，再计算：

\[
v=k_{signed}\frac{\lambda}{2N T_r}.
\]

原始编码谱的 bin 同时含 TX 编码频移，必须在解码时考虑它。

## 3. 图样匹配和复数通道读出

若未编码的目标位于 k，则样例的四个副本位于：

\[
k+[0,208,144,48]\quad (\mathrm{mod}\ 256).
\]

程序逐个扫描候选 k，读取这四处的复数值，并去除各 TX 初始相位。以四处功率的最小值作为图样一致性得分，再做相对阈值和邻域峰值抑制。在选出的候选点上读取四路 TX 复数通道。

图样评分和检测不使用目标真值或真实目标数量。仅在结果验证中计算未编码的参考信号。

按 TX 编码使用的 `np.roll` 是解码的频移补偿；它不改变原始二维 FFT 图的 Bin 顺序。单独移谱并不能隔离完整 TX 时域序列，其他 TX 的副本仍然存在。

## 4. 阅读五张图

| 图 | 内容 |
|---|---|
| 01_ddma_phase_codes.png | 虚构 TX 相位序列及频移 |
| 02_raw_ddma_range_doppler.png | 原始 2D FFT，每个目标有四个编码副本 |
| 03_pattern_alignment.png | 对不同 TX 编码位置进行候选对齐 |
| 04_ddma_decoded_map.png | 图样一致性检测得分，坐标仍为 Bin |
| 05_virtual_channel_validation.png | 复数通道读数与演示参考比较 |

得分图是非线性检测量，不是相干合成后的复数 RD 谱。同距离多个目标可能发生副本碰撞，弱通道也会限制最小功率评分的效果。

默认示例读数如下：

| 初始距离 / 真实速度 | Range Bin | Doppler Bin | 速度读数 |
|---|---:|---:|---:|
| 20 m / −12 m/s | 25 | 194 | −12.08 m/s |
| 45 m / +18 m/s | 58 | 92 | +17.92 m/s |
| 70 m / +70 m/s | 91 | 103 | +20.06 m/s |

第三个目标仍存在原始 PRF 周期的速度混叠。图样匹配帮助辨认 TX 归属，不额外提供解开原始 PRF 混叠的信息。

## 5. 实验

```powershell
# 静止目标：观察编码位置的四个峰
py -3.11 simulate_ddma.py --target 30,0,1 --output ddma_stationary

# 单个运动目标
py -3.11 simulate_ddma.py --target 30,18,1 --output ddma_single

# 实际 chirp 数减半，FFT 仍为 256 点
py -3.11 simulate_ddma.py --chirps 128 --output ddma_short
```

实际 chirp 数和 Doppler FFT 长度需满足代码周期的整倍数要求。实际天线阵列、DOA、射频滤波、定点运算和硬件 CFAR 不在该教学模型中。

## 6. 本地真实配置

真实配置必须放在公开仓库之外，通过 `--config` 指定。脚本自动把其输出放到外部配置旁边的 `local_results/ddma_results`，并拒绝显式写入本公开仓库的输出路径。

公开样例、报告与配图仅使用虚构值。不要将真实配置、真实参数说明或它们生成的配图加入版本控制。

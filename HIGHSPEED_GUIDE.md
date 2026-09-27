# FMCW 教学仿真：完全虚构的公开样例

入口是 [simulate_highspeed_wave.py](simulate_highspeed_wave.py)，公开输入为 [radar_cfg_example.h](configs/radar_cfg_example.h)。本页所有配置数值与结果均来自虚构样例，不对应实际设备。

## 运行

```powershell
py -3.11 simulate_highspeed_wave.py
```

生成 [highspeed_results/REPORT.md](highspeed_results/REPORT.md)、五张图、参数 JSON 和峰值 CSV。

## 样例参数

| 量 | 虚构值 |
|---|---:|
| 起始频率 | 60 GHz |
| 完整上扫带宽 | 300 MHz |
| 上扫 / 回落 / 空闲 | 40 / 4 / 6 μs |
| chirp 起始间隔 | 50 μs |
| ADC 采样率 | 10 MS/s |
| ADC 延迟 / 点数 | 3 μs / 256 点 |
| 距离 FFT / 保留单元 | 256 点 / 前 128 个 |
| 默认实际 chirp 数 / Doppler FFT | 256 / 256 |
| 两维窗口 | Blackman |

调频斜率为 `S = 300 MHz / 40 μs = 7.5 MHz/μs`。ADC 观测时间为 `256 / 10 MHz = 25.6 μs`，因此观测到的扫频带宽为 **192 MHz**。

\[
\Delta R=\frac{c}{2B_{obs}}\approx0.7807\ \mathrm m.
\]

默认距离格点间距也为 0.7807 m。Blackman 窗会拓宽主瓣，因此这个数值不是任意双目标都可分辨的保证。最后一个保留格点约为 99.15 m。

以 ADC 窗口中心频率约 60.1185 GHz 对应的波长计算，默认相干时间为 **12.8 ms**，名义速度分辨尺度约为 **0.1948 m/s**，单通道有符号无模糊速度区间约为 **±24.93 m/s**。

## 回波与处理流程

脚本直接生成解调后的等效复数基带，不对射频载波做离散采样。设目标 `R(t)=R0+v*t`，`tau(t)=2R(t)/c`，采用 `tx * conj(rx)` 混频约定：

\[
\phi(t)=2\pi\left[f_c\tau(t)+St\tau(t)-S\tau(t)^2/2\right]+\phi_0.
\]

正速度表示远离，R0 是第一个 ADC 样本时的距离。模型保留慢时间相位、距离迁移和距离—速度耦合；目标幅度是相对值，另外加入复高斯噪声。

```text
data[256 chirps, 256 ADC samples]
  → ADC 方向加窗和 FFT
range_cube[256 chirps, 128 range bins]
  → chirp 方向加窗和 FFT
range_doppler[256 Doppler bins, 128 range bins]
```

距离 FFT 的复数相位被保留，用于慢时间测速。

## 按顺序看图

1. `01_timing_and_band.png`：扫频、ADC 采样区间、距离频点裁剪。
2. `02_adc_and_range.png`：I/Q 数据及两种窗口下的距离谱。
3. `03_range_doppler.png`：三个演示目标的距离—速度图。
4. `04_slow_time_and_doppler.png`：某距离单元的相位变化及 Doppler 谱。
5. `05_window_and_resolution.png`：格点间距和加窗后的主瓣宽度。

演示目标位于 20、45、70 m，速度分别为 −12、18、70 m/s。最后一个目标会折叠到约 +20.1 m/s。峰值读数采用已知场景数量和邻域抑制，用于教学，不是 CFAR。

## 改参数做实验

```powershell
# 只看一个目标
py -3.11 simulate_highspeed_wave.py --target 30,18,1 --output experiment_single

# 实际观测减半，FFT 仍保留 256 点
py -3.11 simulate_highspeed_wave.py --chirps 128 --output experiment_short

# 实际观测与 FFT 一起翻倍
py -3.11 simulate_highspeed_wave.py --chirps 512 --doppler-fft 512 --output experiment_long
```

补零使频谱格点变密，实际观测时间决定名义速度分辨尺度。模型使用浮点 Blackman 窗，不涉及实际射频滤波器、定点运算、硬件调度或 DOA。

## 使用仅存本地的配置

真实配置应保存在本仓库之外，绝不能覆盖公开样例。用 `--config` 指向外部文件：

```powershell
py -3.11 simulate_highspeed_wave.py --config C:\local-radar\device.local.h
```

外部配置的默认结果位于它旁边的 `local_results/highspeed_results`。即使显式传入 `--output`，脚本也拒绝把外部配置的结果写进公开仓库。头文件内容不被作为 C 代码执行。

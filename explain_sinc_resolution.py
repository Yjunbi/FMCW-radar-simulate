"""用三张图解释：矩形观测窗为什么产生 sinc 谱，以及 1/T 的含义。

运行：py -3.11 explain_sinc_resolution.py
依赖：numpy、matplotlib
输出：plots/01_window.png、02_sinc.png、03_duration.png、04_two_targets.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无需打开图形窗口，也能稳定保存 PNG
import matplotlib.pyplot as plt
import numpy as np


FS = 4000.0  # Hz，采样率
F0 = 40.0  # Hz，单频信号频率
T_SHORT = 0.1  # s，观测时间
T_LONG = 0.2  # s
N_FFT = 80_000  # 仅用于把频谱曲线画平滑；不改变真实分辨率
OUT = Path(__file__).resolve().parent / "plots"


def spectrum(duration: float) -> tuple[np.ndarray, np.ndarray]:
    """计算有限时长复指数信号的归一化频谱。"""
    count = round(FS * duration)
    time = np.arange(count) / FS
    samples = np.exp(2j * np.pi * F0 * time)
    frequency = np.fft.fftshift(np.fft.fftfreq(N_FFT, d=1 / FS))
    magnitude = np.abs(np.fft.fftshift(np.fft.fft(samples, n=N_FFT))) / count
    return frequency, magnitude


def plot_window() -> None:
    """时域：原始单频 × 矩形窗 = 有限时长的观测信号。"""
    time = np.arange(-0.05, 0.2, 1 / FS)
    original = np.cos(2 * np.pi * F0 * time)
    window = ((time >= 0) & (time < T_SHORT)).astype(float)
    observed = original * window

    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    for ax in axes:
        ax.axvspan(0, T_SHORT, color="#d7e8ff", alpha=0.6)
        ax.axvline(0, color="#777777", linestyle=":")
        ax.axvline(T_SHORT, color="#777777", linestyle=":")
        ax.grid(alpha=0.25)
    axes[0].plot(time, original, color="#555555", lw=1.6)
    axes[0].set_ylabel("signal")
    axes[0].set_title("1. Original single tone: cos(2π f₀ t)")
    axes[1].step(time, window, where="post", color="#e28b28", lw=2)
    axes[1].set_ylim(-0.15, 1.25)
    axes[1].set_ylabel("window")
    axes[1].set_title("2. Rectangular window: keep only 0 ≤ t < T")
    axes[2].plot(time, observed, color="#1768ac", lw=1.7)
    axes[2].set_ylabel("observed")
    axes[2].set_xlabel("time (s)")
    axes[2].set_title("3. Observed signal = original × window")
    fig.suptitle(f"Time-domain truncation: f₀ = {F0:g} Hz, T = {T_SHORT:g} s")
    fig.tight_layout()
    fig.savefig(OUT / "01_window.png", dpi=180)
    plt.close(fig)


def plot_sinc() -> None:
    """频域：用计算所得 FFT 验证 sinc 曲线和第一个零点。"""
    frequency, magnitude = spectrum(T_SHORT)
    view = (frequency >= F0 - 4 / T_SHORT) & (frequency <= F0 + 4 / T_SHORT)
    f_view = frequency[view]
    # 连续时间的理论包络：|sin(π(f-f₀)T) / (π(f-f₀)T)|
    ideal = np.abs(np.sinc((f_view - F0) * T_SHORT))

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(f_view, magnitude[view], color="#1768ac", lw=2.2, label="FFT of observed samples")
    ax.plot(f_view, ideal, color="#e28b28", lw=1.5, linestyle="--", label="|sinc((f − f₀)T)|")
    ax.axvline(F0, color="#555555", linestyle=":", label=f"center f₀ = {F0:g} Hz")
    for first_null in (F0 - 1 / T_SHORT, F0 + 1 / T_SHORT):
        ax.axvline(first_null, color="#c33d3d", linestyle="--", alpha=0.85)
        ax.plot(first_null, 0, "o", color="#c33d3d", ms=6)
    ax.annotate("", xy=(F0 + 1 / T_SHORT, 0.72), xytext=(F0, 0.72),
                arrowprops={"arrowstyle": "<->", "color": "#c33d3d", "lw": 1.7})
    ax.text(F0 + 0.5 / T_SHORT, 0.76, "1/T = 10 Hz", ha="center", color="#c33d3d")
    ax.set(xlabel="frequency (Hz)", ylabel="normalized magnitude", ylim=(-0.04, 1.08),
           title="A finite single tone has a sinc-shaped spectrum")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "02_sinc.png", dpi=180)
    plt.close(fig)


def plot_duration() -> None:
    """比较两个观测时间：时间越长，谱峰越窄。"""
    fig, ax = plt.subplots(figsize=(10, 5))
    for duration, color in ((T_SHORT, "#1768ac"), (T_LONG, "#e28b28")):
        frequency, magnitude = spectrum(duration)
        offset = frequency - F0
        view = np.abs(offset) <= 30
        ax.plot(offset[view], magnitude[view], color=color, lw=2,
                label=f"T = {duration:g} s; first null at ±{1 / duration:g} Hz")
        ax.axvline(1 / duration, color=color, linestyle=":", alpha=0.8)
        ax.axvline(-1 / duration, color=color, linestyle=":", alpha=0.8)
    ax.set(xlabel="frequency offset f − f₀ (Hz)", ylabel="normalized magnitude",
           ylim=(-0.04, 1.08), title="Double the observation time → half the main-lobe width")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "03_duration.png", dpi=180)
    plt.close(fig)


def plot_two_targets() -> None:
    """比较两个等强目标的平均功率谱，避免单次相干叠加的相位干涉。"""
    offset = np.linspace(-25, 25, 2001)  # 相对于两目标中点的频率，Hz
    cases = (0.5, 1.0, 1.5)  # 目标间隔 / (1/T)
    fig, axes = plt.subplots(len(cases), 1, figsize=(10, 9), sharex=True, sharey=True)

    for ax, ratio in zip(axes, cases):
        separation = ratio / T_SHORT
        f1_offset, f2_offset = -separation / 2, separation / 2
        power_1 = np.sinc((offset - f1_offset) * T_SHORT) ** 2
        power_2 = np.sinc((offset - f2_offset) * T_SHORT) ** 2
        mean_power = power_1 + power_2

        ax.plot(offset, power_1, color="#1768ac", linestyle="--", lw=1.5,
                label="target 1 alone")
        ax.plot(offset, power_2, color="#e28b28", linestyle="--", lw=1.5,
                label="target 2 alone")
        ax.plot(offset, mean_power, color="#222222", lw=2.4,
                label="sum of powers")
        ax.axvline(f1_offset, color="#1768ac", linestyle=":", alpha=0.8)
        ax.axvline(f2_offset, color="#e28b28", linestyle=":", alpha=0.8)
        ax.set_title(f"Target spacing Δf = {separation:g} Hz = {ratio:g}/T")
        ax.set_ylabel("normalized power")
        ax.set_ylim(-0.05, 1.8)
        ax.grid(alpha=0.25)

    axes[0].legend(loc="upper right", ncol=3, fontsize=9)
    axes[-1].set_xlabel("frequency offset from target midpoint (Hz)")
    fig.suptitle("Two equal-strength targets: phase-averaged power spectra")
    fig.tight_layout()
    fig.savefig(OUT / "04_two_targets.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    plot_window()
    plot_sinc()
    plot_duration()
    plot_two_targets()

    # 数值核对：矩形窗的第一个零点应该出现在 f₀ ± 1/T。
    for duration in (T_SHORT, T_LONG):
        frequency, magnitude = spectrum(duration)
        first_null = F0 + 1 / duration
        index = np.argmin(np.abs(frequency - first_null))
        assert abs(frequency[index] - first_null) < 1e-9
        assert magnitude[index] < 1e-10
        print(f"T={duration:g} s: first null at {first_null:g} Hz; "
              f"normalized magnitude={magnitude[index]:.2e}")
    print(f"Plots saved to: {OUT}")


if __name__ == "__main__":
    main()

"""用五张图学习奈奎斯特采样定理，以及 FMCW 拍频的混叠。

运行：py -3.11 explain_nyquist.py
依赖：numpy、matplotlib
输出：nyquist_plots/ 下的 01～05 PNG 文件。

这里讨论从直流到最高频率 F_max 的实值低通信号。
对这类信号，理想无混叠采样要求 F_s > 2 F_max。
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parent / "nyquist_plots"
BLUE = "#1768ac"
ORANGE = "#e28b28"
RED = "#c43d3d"
GRAY = "#333333"


def save(fig: plt.Figure, filename: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


def plot_enough_samples() -> None:
    """采样足够快时，理想 sinc 插值能恢复带限单频信号。"""
    frequency = 3.0  # Hz
    sample_rate = 20.0  # Hz，远大于 2 * frequency
    t_dense = np.linspace(0, 1, 2001)
    original = np.cos(2 * np.pi * frequency * t_dense)

    # 在显示区间两边也取样，减小有限项 sinc 插值在图边缘的截断误差。
    t_samples_all = np.arange(-2, 3 + 0.5 / sample_rate, 1 / sample_rate)
    samples_all = np.cos(2 * np.pi * frequency * t_samples_all)
    reconstructed = np.sinc(
        sample_rate * (t_dense[:, None] - t_samples_all[None, :])
    ) @ samples_all
    visible = (t_samples_all >= 0) & (t_samples_all <= 1)

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.plot(t_dense, original, color=BLUE, lw=2.3, label="original: 3 Hz cosine")
    ax.plot(t_dense, reconstructed, color=ORANGE, lw=1.7, linestyle="--",
            label="sinc reconstruction")
    ax.scatter(t_samples_all[visible], samples_all[visible], color=RED, s=35,
               zorder=4, label="samples: Fs = 20 Hz")
    ax.set(xlabel="time (s)", ylabel="amplitude", ylim=(-1.2, 1.2),
           title="1. Enough samples: the waveform can be reconstructed")
    ax.grid(alpha=0.25)
    ax.legend(ncol=3, loc="upper right")
    save(fig, "01_enough_samples.png")
    print(f"Figure 1: maximum reconstruction error in view = "
          f"{np.max(np.abs(reconstructed - original)):.3g}")


def plot_aliasing() -> None:
    """3 Hz 与 7 Hz 余弦在 10 Hz 采样下产生完全相同的样本。"""
    sample_rate = 10.0
    t_dense = np.linspace(0, 1, 3001)
    t_samples = np.arange(11) / sample_rate
    low = np.cos(2 * np.pi * 3 * t_dense)
    high = np.cos(2 * np.pi * 7 * t_dense)
    samples_3 = np.cos(2 * np.pi * 3 * t_samples)
    samples_7 = np.cos(2 * np.pi * 7 * t_samples)

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.plot(t_dense, low, color=BLUE, lw=2.2, label="3 Hz cosine")
    ax.plot(t_dense, high, color=ORANGE, lw=1.7, linestyle="--", label="7 Hz cosine")
    ax.scatter(t_samples, samples_3, color=RED, edgecolor="white", s=55,
               zorder=4, label="shared samples at Fs = 10 Hz")
    ax.set(xlabel="time (s)", ylabel="amplitude", ylim=(-1.2, 1.2),
           title="2. Aliasing: two different waves give identical samples")
    ax.grid(alpha=0.25)
    ax.legend(ncol=3, loc="upper right")
    save(fig, "02_aliasing.png")
    print(f"Figure 2: largest sample difference (3 Hz vs 7 Hz) = "
          f"{np.max(np.abs(samples_3 - samples_7)):.2e}")


def triangular_spectrum(frequency: np.ndarray, center: float, f_max: float) -> np.ndarray:
    """示意性带限频谱；只在 [center-F_max, center+F_max] 非零。"""
    return np.maximum(1 - np.abs(frequency - center) / f_max, 0)


def plot_spectral_copies() -> None:
    """采样后频谱每隔 Fs 复制一次；副本重叠即混叠。"""
    f_max = 3.0
    frequency = np.linspace(-12, 12, 2401)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True, sharey=True)

    for ax, sample_rate in zip(axes, (8.0, 5.0)):
        base = triangular_spectrum(frequency, 0, f_max)
        left = triangular_spectrum(frequency, -sample_rate, f_max)
        right = triangular_spectrum(frequency, sample_rate, f_max)
        ax.fill_between(frequency, base, color=BLUE, alpha=0.20)
        ax.fill_between(frequency, left, color=ORANGE, alpha=0.20)
        ax.fill_between(frequency, right, color=ORANGE, alpha=0.20)
        ax.plot(frequency, base, color=BLUE, lw=2, label="original spectrum")
        ax.plot(frequency, left, color=ORANGE, lw=1.7, linestyle="--",
                label="shifted copies")
        ax.plot(frequency, right, color=ORANGE, lw=1.7, linestyle="--")
        ax.plot(frequency, base + left + right, color=GRAY, lw=1.5,
                label="sampled spectrum (sum)")
        ax.set_ylabel("relative spectrum")
        ax.set_title(f"Fs = {sample_rate:g} Hz; 2 Fmax = {2 * f_max:g} Hz: "
                     + ("no overlap" if sample_rate > 2 * f_max else "overlap = aliasing"))
        ax.set_ylim(0, 1.15)
        ax.grid(alpha=0.2)

    # Fs=5、Fmax=3 时，中央频谱和相邻副本在 [2,3] 与 [-3,-2] 重叠。
    for left_edge, right_edge in ((-3, -2), (2, 3)):
        axes[1].axvspan(left_edge, right_edge, color=RED, alpha=0.13)
    axes[0].legend(ncol=3, loc="upper right")
    axes[1].set_xlabel("frequency (Hz)")
    fig.suptitle("3. Sampling creates spectral copies spaced by Fs")
    save(fig, "03_spectral_copies.png")


def plot_boundary() -> None:
    """在 Fs=2f 的边界，正弦采样点可以全部为零。"""
    frequency = 5.0
    sample_rate = 10.0
    t_dense = np.linspace(0, 0.6, 2001)
    t_samples = np.arange(7) / sample_rate
    original = np.sin(2 * np.pi * frequency * t_dense)
    samples = np.sin(2 * np.pi * frequency * t_samples)

    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.plot(t_dense, original, color=BLUE, lw=2, label="5 Hz sine")
    ax.axhline(0, color=GRAY, linestyle=":", label="zero signal")
    ax.scatter(t_samples, samples, color=RED, s=60, zorder=4,
               label="samples at Fs = 10 Hz: all zero")
    ax.set(xlabel="time (s)", ylabel="amplitude", ylim=(-1.2, 1.2),
           title="4. At the exact boundary Fs = 2f, phase can make samples ambiguous")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    save(fig, "04_boundary.png")
    print(f"Figure 4: largest boundary-sample magnitude = {np.max(np.abs(samples)):.2e}")


def plot_radar_beat_alias() -> None:
    """实值 ADC 的 FMCW 拍频超过 Fs/2 后，会折回低频。"""
    sample_rate = 10.0  # MHz
    true_beat = np.linspace(0, 25, 2501)  # MHz
    # 实值余弦采样后，把离散频率折叠到 [0, Fs/2]。
    observed = np.abs((true_beat + sample_rate / 2) % sample_rate - sample_rate / 2)

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.axvspan(sample_rate / 2, true_beat[-1], color=RED, alpha=0.07,
               label="above real-ADC Nyquist frequency")
    ax.plot(true_beat, observed, color=BLUE, lw=2.2, label="observed positive frequency")
    ax.plot([0, 5], [0, 5], color=GRAY, linestyle="--", lw=1.2,
            label="unaliased region")
    ax.scatter([7], [3], color=RED, s=60, zorder=4)
    ax.annotate("7 MHz → 3 MHz", xy=(7, 3), xytext=(9.3, 4.1),
                arrowprops={"arrowstyle": "->", "color": RED}, color=RED)
    ax.axvline(sample_rate / 2, color=RED, linestyle=":")
    ax.set(xlabel="true beat frequency (MHz)", ylabel="observed frequency (MHz)",
           xlim=(0, 25), ylim=(-0.2, 5.5),
           title="5. FMCW example: a real ADC folds beat frequencies above Fs/2")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    save(fig, "05_radar_beat_alias.png")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    plot_enough_samples()
    plot_aliasing()
    plot_spectral_copies()
    plot_boundary()
    plot_radar_beat_alias()
    print(f"Five figures saved to: {OUT}")


if __name__ == "__main__":
    main()

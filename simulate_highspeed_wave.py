"""使用完全虚构的 radar_cfg_example.h 学习 FMCW 距离/速度处理。

运行：py -3.11 simulate_highspeed_wave.py
改场景：py -3.11 simulate_highspeed_wave.py --target 30,20,1 --target 60,-15,0.7
改实际 chirp 数：py -3.11 simulate_highspeed_wave.py --chirps 128

模型约定：MHz/微秒单位；等间隔 chirp；正速度表示远离雷达；
等效复数单通道。默认选取与示例 Doppler FFT 等长的实际 chirp 序列，
这是教学场景的选择。外部配置的输出自动放在仓库外。程序不执行 C 代码。
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C = 299_792_458.0
ROOT = Path(__file__).resolve().parent
PUBLIC_CONFIG = ROOT / "configs" / "radar_cfg_example.h"


def choose_output(config: Path, requested: Path | None, name: str) -> Path:
    """公共样例可生成仓库内配图；其他配置的结果必须存到仓库外。"""
    example = config.resolve() == PUBLIC_CONFIG.resolve()
    output = requested if requested is not None else (
        ROOT / name if example else config.resolve().parent / "local_results" / name
    )
    output = output.resolve()
    if not example and (output == ROOT or ROOT in output.parents):
        raise ValueError("Custom configurations require an output directory outside this public repository.")
    return output


def field(text: str, name: str) -> str:
    match = re.search(r"\." + re.escape(name) + r"\s*=\s*([^,\n\\]+)", text)
    if not match:
        raise ValueError(f"Header field not found: {name}")
    return match.group(1).strip()


def number(text: str, name: str) -> float:
    value = re.sub(r"^\(u32\)\s*", "", field(text, name))
    value = re.sub(r"[uUfF]$", "", value)
    return float(value)


@dataclass(frozen=True)
class Config:
    start_hz: float
    bandwidth_hz: float
    rise_s: float
    fall_s: float
    idle_s: float
    fs_hz: float
    samples: int
    adc_delay_s: float
    range_fft: int
    range_bins: int
    doppler_fft: int

    @classmethod
    def read(cls, path: Path) -> Config:
        text = path.read_text(encoding="utf-8-sig")
        adc = re.fullmatch(r"RADAR_ADC_RATE_(\d+)M", field(text, "adc_sample_rate"))
        if adc is None:
            raise ValueError("Unsupported ADC rate enum; update its conversion explicitly.")
        for dimension in ("range", "doppler"):
            if number(text, f"{dimension}_win_en") != 1 or field(text, f"{dimension}_win_type") != "WIN_BLACKMAN":
                raise ValueError("This teaching script expects enabled Blackman windows.")
        cfg = cls(
            number(text, "chirp_startfreq_khz") * 1e3,
            number(text, "chirp_bandwidth") * 1e6,
            number(text, "chirp_rise") * 1e-6,
            number(text, "chirp_fall") * 1e-6,
            number(text, "chirp_idle") * 1e-6,
            float(adc.group(1)) * 1e6,
            int(number(text, "adc_chirp_size")),
            number(text, "adc_sample_delay") * 1e-6,
            int(number(text, "range_fft_size")),
            int(number(text, "num_rangebin")),
            int(number(text, "doppler_fft_size")),
        )
        if cfg.adc_delay_s + cfg.adc_time > cfg.rise_s + 1e-12:
            raise ValueError("Assumed ADC window extends past the rising ramp.")
        if cfg.range_fft < cfg.samples or not 0 < cfg.range_bins <= cfg.range_fft // 2:
            raise ValueError("Unsupported range FFT/crop dimensions.")
        return cfg

    @property
    def slope(self) -> float:
        return self.bandwidth_hz / self.rise_s

    @property
    def pri(self) -> float:
        return self.rise_s + self.fall_s + self.idle_s

    @property
    def adc_time(self) -> float:
        return self.samples / self.fs_hz

    @property
    def sampled_bandwidth(self) -> float:
        return self.slope * self.adc_time

    @property
    def reference_hz(self) -> float:
        # 采用 ADC 窗口中心的 RF 频率作为速度坐标的参考频率。
        return self.start_hz + self.slope * (self.adc_delay_s + self.adc_time / 2)

    @property
    def wavelength(self) -> float:
        return C / self.reference_hz

    @property
    def range_step(self) -> float:
        return C * self.fs_hz / (2 * self.slope * self.range_fft)


@dataclass(frozen=True)
class Target:
    range_m: float
    velocity_m_s: float
    amplitude: float = 1.0
    phase_rad: float = 0.0


def target_arg(value: str) -> Target:
    try:
        parts = [float(x) for x in value.split(",")]
        if len(parts) not in (2, 3):
            raise ValueError
        return Target(*parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use range_m,velocity_m_s[,amplitude].") from exc


def alias_velocity(cfg: Config, velocity: float) -> float:
    period = cfg.wavelength / (2 * cfg.pri)
    return (velocity + period / 2) % period - period / 2


def prediction(cfg: Config, target: Target, chirps: int) -> dict:
    t_fast = cfg.adc_delay_s + (cfg.samples - 1) / (2 * cfg.fs_hz)
    t_center = (chirps - 1) * cfg.pri / 2 + (cfg.samples - 1) / (2 * cfg.fs_hz)
    range_center = target.range_m + target.velocity_m_s * t_center
    tau = 2 * range_center / C
    # 与生成回波同一模型的相位导数；包括 Doppler 对快时间拍频的影响。
    doppler_hz = (cfg.start_hz + cfg.slope * t_fast - cfg.slope * tau) * 2 * target.velocity_m_s / C
    beat_hz = cfg.slope * tau + doppler_hz
    doppler_alias_hz = (doppler_hz + 1 / (2 * cfg.pri)) % (1 / cfg.pri) - 1 / (2 * cfg.pri)
    return {
        "range_at_first_adc_m": target.range_m,
        "range_at_cpi_center_m": range_center,
        "physical_velocity_m_s": target.velocity_m_s,
        "predicted_apparent_range_m": C * beat_hz / (2 * cfg.slope),
        "predicted_aliased_velocity_m_s": doppler_alias_hz * cfg.wavelength / 2,
        "range_migration_m": target.velocity_m_s * (chirps - 1) * cfg.pri,
    }


def simulate(cfg: Config, targets: list[Target], chirps: int, snr_db: float | None,
             seed: int) -> np.ndarray:
    """生成 [chirp, ADC sample] 数据。R0 对应第一个 ADC 样本的时刻。

    tx(t) = exp(j2π(fc*t + S*t²/2))；rx 是延迟 tau 的副本。
    tx * conj(rx) 的相位 = 2π(fc*tau + S*t*tau - S*tau²/2)。
    tau 随目标运动变化，因而同时保留慢时间 Doppler、距离迁移和拍频耦合。
    使用 R(t)/c 的常见非相对论延迟近似；假定各 chirp 相干，忽略传播损耗。
    """
    fast = cfg.adc_delay_s + np.arange(cfg.samples) / cfg.fs_hz
    elapsed = np.arange(chirps)[:, None] * cfg.pri + np.arange(cfg.samples)[None, :] / cfg.fs_hz
    data = np.zeros((chirps, cfg.samples), dtype=np.complex128)
    for target in targets:
        distance = target.range_m + target.velocity_m_s * elapsed
        tau = 2 * distance / C
        if np.min(distance) <= 0 or np.max(tau) >= cfg.adc_delay_s:
            raise ValueError("Target must remain positive and its echo must stay within the same ramp.")
        cycles = (cfg.start_hz + cfg.slope * fast[None, :]) * tau - 0.5 * cfg.slope * tau**2
        data += target.amplitude * np.exp(1j * (2 * np.pi * cycles + target.phase_rad))
    if snr_db is not None:
        rng = np.random.default_rng(seed)
        variance = sum(t.amplitude**2 for t in targets) / 10**(snr_db / 10)
        data += np.sqrt(variance / 2) * (rng.standard_normal(data.shape) + 1j * rng.standard_normal(data.shape))
    return data


def process(cfg: Config, data: np.ndarray, doppler_fft: int) -> tuple:
    if doppler_fft < len(data):
        raise ValueError("Doppler FFT must be >= actual chirps; truncation would discard observations.")
    range_window = np.blackman(cfg.samples)
    slow_window = np.blackman(len(data))
    range_cube = np.fft.fft(data * range_window, n=cfg.range_fft, axis=1)[:, :cfg.range_bins] / range_window.sum()
    rd = np.fft.fftshift(np.fft.fft(range_cube * slow_window[:, None], n=doppler_fft, axis=0), axes=0) / slow_window.sum()
    ranges = np.arange(cfg.range_bins) * cfg.range_step
    velocities = np.fft.fftshift(np.fft.fftfreq(doppler_fft, cfg.pri)) * cfg.wavelength / 2
    return range_cube, rd, ranges, velocities


def db(value: np.ndarray, reference: float | None = None) -> np.ndarray:
    magnitude = np.abs(value)
    if reference is None:
        reference = float(magnitude.max())
    return 20 * np.log10(np.maximum(magnitude / max(reference, 1e-30), 1e-6))


def read_peaks(rd: np.ndarray, ranges: np.ndarray, velocities: np.ndarray, count: int) -> list[dict]:
    """仅用于教学的全局峰值选取和邻域抑制，不声称复刻配置中的 CFAR。"""
    remaining = np.abs(rd).copy()
    peaks = []
    for _ in range(count):
        d, r = np.unravel_index(np.argmax(remaining), remaining.shape)
        peaks.append({"range_m": float(ranges[r]), "velocity_m_s": float(velocities[d]),
                      "magnitude": float(remaining[d, r])})
        for row in range(d - 4, d + 5):
            remaining[row % len(velocities), max(0, r - 4):r + 5] = 0
    return peaks


def save(fig: plt.Figure, output: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(output / name, dpi=160)
    plt.close(fig)


def draw(cfg: Config, targets: list[Target], data: np.ndarray, cube: np.ndarray,
         rd: np.ndarray, ranges: np.ndarray, velocities: np.ndarray,
         predictions: list[dict], peaks: list[dict], output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    for m in range(3):
        start = m * cfg.pri
        axes[0].plot(np.array([start, start + cfg.rise_s]) * 1e6,
                     [0, cfg.bandwidth_hz / 1e6], color="#1768ac", lw=2)
        axes[0].plot(np.array([start + cfg.rise_s, start + cfg.rise_s + cfg.fall_s, start + cfg.pri]) * 1e6,
                     [cfg.bandwidth_hz / 1e6, 0, 0], color="#777777", linestyle="--")
        axes[0].axvspan((start + cfg.adc_delay_s) * 1e6,
                       (start + cfg.adc_delay_s + cfg.adc_time) * 1e6,
                       color="#e28b28", alpha=0.22, label="ADC observation" if m == 0 else None)
    axes[0].set(xlabel="time (us)", ylabel=f"RF offset from {cfg.start_hz / 1e9:g} GHz (MHz)",
                title=f"Assumed PRI = rise + fall + idle = {cfg.pri * 1e6:.3f} us")
    axes[0].legend()
    axes[1].axvspan(0, cfg.range_bins * cfg.fs_hz / cfg.range_fft / 1e6, color="#1768ac", alpha=0.25,
                   label=f"retained example bins 0..{cfg.range_bins - 1}")
    axes[1].axvline(cfg.fs_hz / 2e6, color="#c43d3d", linestyle=":", label="Fs/2 (real-ADC reference)")
    axes[1].set(xlim=(0, cfg.fs_hz / 2e6 + 2), ylim=(0, 1), yticks=[], xlabel="positive beat frequency (MHz)",
                title="Synthetic example: FFT crop and real-ADC Nyquist reference")
    axes[1].legend(loc="upper left")
    for ax in axes:
        ax.grid(alpha=0.2)
    save(fig, output, "01_timing_and_band.png")

    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    shown = min(100, cfg.samples)
    time_us = np.arange(shown) / cfg.fs_hz * 1e6
    axes[0].plot(time_us, data[0, :shown].real, label="I = real part")
    axes[0].plot(time_us, data[0, :shown].imag, label="Q = imaginary part", alpha=0.8)
    axes[0].set(xlabel="time since ADC start (us)", ylabel="normalized amplitude",
                title="First chirp: a mixture of target beat signals + noise")
    axes[0].legend()
    rectangular = np.fft.fft(data[0], n=cfg.range_fft)[:cfg.range_bins] / cfg.samples
    reference = max(np.max(np.abs(rectangular)), np.max(np.abs(cube[0])))
    axes[1].plot(ranges, db(rectangular, reference), label="rectangular window", alpha=0.8)
    axes[1].plot(ranges, db(cube[0], reference), label="configured Blackman (float approximation)")
    axes[1].set(xlabel="apparent range (m)", ylabel="relative magnitude (dB)", ylim=(-80, 5),
                title="Range FFT: Blackman reduces sidelobes and broadens peaks")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=0.25)
    save(fig, output, "02_adc_and_range.png")

    fig, ax = plt.subplots(figsize=(11, 7))
    mesh = ax.pcolormesh(ranges, velocities, db(rd), shading="auto", cmap="magma", vmin=-65, vmax=0)
    fig.colorbar(mesh, ax=ax, label="relative magnitude (dB)")
    for i, (target, pred) in enumerate(zip(targets, predictions), 1):
        r, v = pred["predicted_apparent_range_m"], pred["predicted_aliased_velocity_m_s"]
        ax.plot(r, v, "o", ms=10, mfc="none", mec="cyan", mew=1.5)
        ax.annotate(f"T{i}: physical v={target.velocity_m_s:g} m/s", (r, v), xytext=(7, 10),
                    textcoords="offset points", color="white", fontsize=9)
    ax.scatter([p["range_m"] for p in peaks], [p["velocity_m_s"] for p in peaks], marker="x",
               color="lime", s=40, label="simple peak readout (not hardware CFAR)")
    ax.set(xlabel="apparent range (m)", ylabel="aliased radial velocity (m/s)",
           title=f"Range-Doppler: {len(data)} actual chirps; Blackman in both dimensions")
    ax.legend(loc="upper right", fontsize=9)
    save(fig, output, "03_range_doppler.png")

    # 选速度绝对值最小的目标，用其距离单元解释慢时间相位。
    chosen = min(range(len(targets)), key=lambda i: abs(targets[i].velocity_m_s))
    pred = predictions[chosen]
    r_bin = int(np.argmin(np.abs(ranges - pred["predicted_apparent_range_m"])))
    slow = cube[:, r_bin]
    time_ms = np.arange(len(data)) * cfg.pri * 1e3
    phase = np.unwrap(np.angle(slow))
    expected_fd = pred["predicted_aliased_velocity_m_s"] * 2 / cfg.wavelength
    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    axes[0].plot(time_ms, phase - phase[0], label="unwrapped phase of selected range bin")
    axes[0].plot(time_ms, 2 * np.pi * expected_fd * time_ms / 1e3, "--", label="predicted phase slope")
    axes[0].set(xlabel="slow time (ms)", ylabel="phase change (rad)",
                title=f"Target {chosen + 1}: phase evolution across chirps")
    axes[1].plot(velocities, db(rd[:, r_bin]))
    axes[1].axvline(pred["predicted_aliased_velocity_m_s"], color="#e28b28", linestyle="--", label="predicted Doppler velocity")
    axes[1].set(xlabel="radial velocity (m/s)", ylabel="relative magnitude (dB)", ylim=(-85, 5),
                title=f"Doppler FFT at range bin {r_bin}")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend()
    save(fig, output, "04_slow_time_and_doppler.png")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, count, interval, conversion, label in (
        (axes[0], cfg.samples, 1 / cfg.fs_hz, C / (2 * cfg.slope), "range offset (m)"),
        (axes[1], len(data), cfg.pri, cfg.wavelength / 2, "velocity offset (m/s)"),
    ):
        fine_fft = max(65536, count * 128)
        axis = np.fft.fftshift(np.fft.fftfreq(fine_fft, interval)) * conversion
        nominal = conversion / (count * interval)
        visible = np.abs(axis) <= 6 * nominal
        for window, name in ((np.ones(count), "rectangular"), (np.blackman(count), "Blackman")):
            response = np.fft.fftshift(np.fft.fft(window, n=fine_fft))
            ax.plot(axis[visible], db(response)[visible], label=name)
        ax.axvline(nominal, color="#777777", linestyle=":", label=f"nominal = {nominal:.3f}")
        ax.set(xlabel=label, ylabel="normalized response (dB)", ylim=(-90, 3))
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle("FFT grid spacing is not Blackman two-target resolving ability")
    save(fig, output, "05_window_and_resolution.png")


def verify(cfg: Config, chirps: int, doppler_fft: int) -> dict:
    """用独立场景检查静止目标、运动目标和混叠目标的坐标，容许一个格点。"""
    results = {}
    for name, target in (("stationary", Target(30, 0)), ("moving", Target(30, 15)),
                         ("velocity_alias", Target(30, cfg.wavelength / (4 * cfg.pri) + 8))):
        data = simulate(cfg, [target], chirps, None, 0)
        _, rd, ranges, velocities = process(cfg, data, doppler_fft)
        d, r = np.unravel_index(np.argmax(np.abs(rd)), rd.shape)
        expected = prediction(cfg, target, chirps)
        range_error = abs(ranges[r] - expected["predicted_apparent_range_m"])
        period = cfg.wavelength / (2 * cfg.pri)
        speed_error = abs((velocities[d] - expected["predicted_aliased_velocity_m_s"] + period / 2) % period - period / 2)
        assert range_error <= cfg.range_step, (name, range_error)
        assert speed_error <= period / doppler_fft, (name, speed_error)
        results[name] = {"range_error_m": float(range_error), "velocity_error_m_s": float(speed_error)}
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=PUBLIC_CONFIG)
    parser.add_argument("--output", type=Path, default=None, help="Custom configurations must write outside the repository.")
    parser.add_argument("--chirps", type=int, default=None, help="Actual coherent samples; the teaching default equals the configured FFT length.")
    parser.add_argument("--doppler-fft", type=int, default=None, help="Override configured FFT length; must be >= chirps.")
    parser.add_argument("--snr-db", type=float, default=15.0, help="Aggregate signal/noise power ratio per complex ADC sample.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--target", type=target_arg, action="append", help="range_m,velocity_m_s[,amplitude]; repeat for multiple targets")
    args = parser.parse_args()
    cfg = Config.read(args.config)
    chirps = cfg.doppler_fft if args.chirps is None else args.chirps
    try:
        args.output = choose_output(args.config, args.output, "highspeed_results")
    except ValueError as exc:
        parser.error(str(exc))
    doppler_fft = cfg.doppler_fft if args.doppler_fft is None else args.doppler_fft
    if not 8 <= chirps <= 8192 or doppler_fft < chirps:
        parser.error("Use 8..8192 actual chirps and a Doppler FFT length >= actual chirps.")
    targets = args.target or [Target(20, -12, 1.0, 0.2), Target(45, 18, 0.8, 1.0), Target(70, 70, 0.6, 2.0)]
    if any(not np.isfinite([t.range_m, t.velocity_m_s, t.amplitude]).all() or t.amplitude <= 0 for t in targets):
        parser.error("Targets must have finite values and positive amplitude.")
    predictions = [prediction(cfg, target, chirps) for target in targets]
    if any(not 0 <= p["predicted_apparent_range_m"] < (cfg.range_bins - 1) * cfg.range_step for p in predictions):
        parser.error("Choose targets whose predicted range peaks stay in the retained FFT bins.")
    args.output.mkdir(parents=True, exist_ok=True)
    data = simulate(cfg, targets, chirps, args.snr_db, args.seed)
    cube, rd, ranges, velocities = process(cfg, data, doppler_fft)
    peaks = read_peaks(rd, ranges, velocities, len(targets))
    checks = verify(cfg, chirps, doppler_fft)
    draw(cfg, targets, data, cube, rd, ranges, velocities, predictions, peaks, args.output)

    metrics = {
        "chirp_interval_us_assumed": cfg.pri * 1e6,
        "slope_MHz_per_us": cfg.slope / 1e12,
        "adc_observation_us": cfg.adc_time * 1e6,
        "observed_bandwidth_MHz": cfg.sampled_bandwidth / 1e6,
        "full_sweep_nominal_range_resolution_m": C / (2 * cfg.bandwidth_hz),
        "observed_nominal_range_resolution_m": C / (2 * cfg.sampled_bandwidth),
        "range_grid_m": cfg.range_step,
        "last_retained_range_bin_m": float(ranges[-1]),
        "reference_RF_GHz": cfg.reference_hz / 1e9,
        "assumed_actual_chirps": chirps,
        "configured_doppler_fft": cfg.doppler_fft,
        "used_doppler_fft": doppler_fft,
        "coherent_time_ms_assumed": chirps * cfg.pri * 1e3,
        "nominal_velocity_resolution_m_s": cfg.wavelength / (2 * chirps * cfg.pri),
        "velocity_grid_m_s": cfg.wavelength / (2 * doppler_fft * cfg.pri),
        "max_signed_velocity_abs_m_s": cfg.wavelength / (4 * cfg.pri),
    }
    assumptions = [
        "chirp_bandwidth interpreted as MHz; rise/fall/idle/ADC delay interpreted as us.",
        "PRI assumed to be rise + fall + idle, with no extra gaps inside this CPI.",
        "Actual chirps explicitly selected by --chirps; FFT length is not evidence of actual chirp count.",
        "A coherent, dechirped complex equivalent channel is simulated. Actual ADC format is unspecified.",
        "The bundled configuration is entirely synthetic and contains no physical device settings.",
        "This single-channel model does not apply the example TX phase codes; use simulate_ddma.py for coding.",
        "The CPI length is selected for the simulation; no hardware frame schedule is inferred.",
        f"num_rangebin interpreted as keeping bins 0..{cfg.range_bins - 1}.",
        "Blackman windows use numpy symmetric floating coefficients, not hardware quantized tables.",
        "RF filters, gains, fixed-point scaling, CFAR/NVE and DOA are not emulated.",
        "Targets have fixed relative amplitudes, constant radial speeds, no acceleration or multipath.",
    ]
    report = {"source_header": args.config.name, "source_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
              "interpreted_config_SI": asdict(cfg), "metrics": metrics, "assumptions": assumptions,
              "targets": [asdict(t) for t in targets], "predictions": predictions, "simple_peaks": peaks,
              "numerical_checks": checks, "snr_db": args.snr_db, "seed": args.seed}
    (args.output / "config_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    with (args.output / "peak_estimates.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(peaks[0]))
        writer.writeheader()
        writer.writerows(peaks)
    lines = ["# 本次仿真结果", "", "常规 FMCW 单通道教学仿真。仓库自带配置使用完全虚构的参数。", "",
             "## 参数换算", "", "| 参数 | 数值 |", "|---|---:|"]
    lines += [f"| {key} | {value:.9g} |" for key, value in metrics.items()]
    lines += ["", "## 场景与预测", "", "| 目标 | 初始距离 m | 真实速度 m/s | 预测表观距离 m | 预测混叠速度 m/s |",
              "|---|---:|---:|---:|---:|"]
    lines += [f"| T{i} | {t.range_m:g} | {t.velocity_m_s:g} | {p['predicted_apparent_range_m']:.3f} | {p['predicted_aliased_velocity_m_s']:.3f} |"
              for i, (t, p) in enumerate(zip(targets, predictions), 1)]
    lines += ["", "表观距离含 Doppler 耦合和 CPI 内运动；速度坐标采用 ADC 中心 RF 波长。",
              "峰值读数只有 FFT 格点精度，采用简单峰值选取；不等同于配置的硬件 CFAR。", "",
              "## 阅读图片", ""]
    for name in ("01_timing_and_band", "02_adc_and_range", "03_range_doppler", "04_slow_time_and_doppler", "05_window_and_resolution"):
        lines += [f"![{name}]({name}.png)", ""]
    lines += ["## 模型约定", ""] + [f"- {x}" for x in assumptions]
    (args.output / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    print("Numerical checks passed: stationary target, moving target, velocity alias.")
    print(f"Plots and report: {args.output.resolve()}")


if __name__ == "__main__":
    main()

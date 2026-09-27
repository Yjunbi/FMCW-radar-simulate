"""4TX/1RX DDMA 教学仿真：相位编码、谱搬移、图样匹配和 TX 复数通道读出。

运行：py -3.11 simulate_ddma.py
单目标：py -3.11 simulate_ddma.py --target 30,18,1

复用 simulate_highspeed_wave.py 的波形和运动模型。默认假设 4TX 同时发射、
256 个相位码对应 2π，并连续累加相位。使用 tx_ref * conj(rx) 混频，
因此发射端正相位步进在基带产生负多普勒搬移；--code-sign 可切换约定。
图样匹配适用于谱副本可以辨认的稀疏场景，不保证解决同距离多目标的所有碰撞。
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from simulate_highspeed_wave import Config, ROOT, PUBLIC_CONFIG, Target, choose_output, db, prediction, save, simulate, target_arg


def read_codes(text: str, name: str) -> np.ndarray:
    match = re.search(r"\." + re.escape(name) + r"\s*=\s*\{([^}]+)\}", text)
    if match is None:
        raise ValueError(f"Cannot read {name} from header.")
    tokens = match.group(1).replace("\\", "").split(",")
    return np.array([int(re.sub(r"[uUlL]+$", "", x.strip()), 0) for x in tokens if x.strip()], dtype=int)


@dataclass
class DDMA:
    starts: np.ndarray
    steps: np.ndarray
    levels: int
    sign: int

    @property
    def period(self) -> int:
        return math.lcm(*(self.levels // math.gcd(int(x), self.levels) for x in self.steps))

    @property
    def offsets(self) -> np.ndarray:
        """折到 [-0.5,0.5) 的有效基带频移，单位 cycles/chirp。"""
        return (self.sign * self.steps / self.levels + 0.5) % 1 - 0.5

    @property
    def initial_phase(self) -> np.ndarray:
        return self.sign * 2 * np.pi * self.starts / self.levels

    def codes(self, chirps: int) -> np.ndarray:
        phase_code = (self.starts[:, None] + self.steps[:, None] * np.arange(chirps)) % self.levels
        return np.exp(1j * self.sign * 2 * np.pi * phase_code / self.levels)

    def shifts(self, fft_size: int) -> np.ndarray:
        exact = self.offsets * fft_size
        if not np.allclose(exact, np.round(exact), atol=1e-10):
            raise ValueError("Doppler FFT size must place all code offsets on integer bins.")
        return np.round(exact).astype(int)

    def rotational_symmetries(self, fft_size: int) -> list[int]:
        occupied = set((self.shifts(fft_size) % fft_size).tolist())
        return [shift for shift in range(fft_size)
                if {(x + shift) % fft_size for x in occupied} == occupied]


def encode(base: np.ndarray, coding: DDMA, gains: np.ndarray, snr_db: float | None,
           seed: int) -> np.ndarray:
    """四个 TX 的已编码回波先相加，接收端只获得这一份混合数据。"""
    codes = coding.codes(len(base))
    mixture = base * np.sum(gains[:, None] * codes, axis=0)[:, None]
    if snr_db is not None:
        rng = np.random.default_rng(seed)
        variance = np.mean(np.abs(mixture)**2) / 10**(snr_db / 10)
        mixture += np.sqrt(variance / 2) * (rng.standard_normal(base.shape) + 1j * rng.standard_normal(base.shape))
    return mixture


def process_ddma(cfg: Config, data: np.ndarray, fft_size: int) -> tuple:
    """保留原始 FFT 输出顺序；二维图直接使用 Range Bin / Doppler Bin。

    Doppler Bin 0 是 DC，正频率在前、负频率在后。仅用 fftfreq 换算读数，
    不对 FFT 数组重新排序。TX 编码补偿在 align_and_score 中单独执行。
    """
    range_window = np.blackman(cfg.samples)
    slow_window = np.blackman(len(data))
    range_cube = np.fft.fft(data * range_window, n=cfg.range_fft, axis=1)[:, :cfg.range_bins] / range_window.sum()
    rd = np.fft.fft(range_cube * slow_window[:, None], n=fft_size, axis=0) / slow_window.sum()
    ranges = np.arange(cfg.range_bins) * cfg.range_step
    velocities = np.fft.fftfreq(fft_size, cfg.pri) * cfg.wavelength / 2
    return range_cube, rd, ranges, velocities


def align_and_score(rd: np.ndarray, coding: DDMA) -> tuple[np.ndarray, np.ndarray]:
    """对每个物理 Doppler 假设取四个编码位置，去除初始编码相位。

    aligned[q,d,r] = Y[d + shift_q,r] * exp(-j*initial_phase_q)。
    移谱本身还没有分离 TX，其他 TX 的副本仍在；取四支路的最小功率作为
    一致性得分，要求四个预测位置都有响应。这个得分不是相干合成的信号。
    """
    aligned = np.stack([np.roll(rd, -int(shift), axis=0) * np.exp(-1j * phase)
                        for shift, phase in zip(coding.shifts(len(rd)), coding.initial_phase)])
    score = np.min(np.abs(aligned)**2, axis=0)
    return aligned, score


def detect(score: np.ndarray, threshold_db: float, max_peaks: int = 20) -> list[tuple[int, int]]:
    """仅用接收图样得分做相对阈值与邻域抑制，不读取目标真值和目标数量。"""
    work = score.copy()
    threshold = float(work.max()) * 10**(threshold_db / 10)
    detections = []
    while len(detections) < max_peaks:
        d, r = np.unravel_index(np.argmax(work), work.shape)
        if work[d, r] < threshold or work[d, r] <= 0:
            break
        detections.append((int(d), int(r)))
        for row in range(d - 4, d + 5):
            work[row % work.shape[0], max(0, r - 4):r + 5] = 0
    return detections


def wrap_velocity(velocity: float | np.ndarray, period: float) -> float | np.ndarray:
    return (velocity + period / 2) % period - period / 2


def numerical_checks(cfg: Config, coding: DDMA, chirps: int, fft_size: int,
                     gains: np.ndarray) -> dict:
    # 验证相位调制与整数格点循环移谱等价，包括非零初始相位。
    rng = np.random.default_rng(19)
    sequence = (rng.standard_normal(chirps) + 1j * rng.standard_normal(chirps)) * np.blackman(chirps)
    base_spectrum = np.fft.fft(sequence, n=fft_size)
    shift_error = 0.0
    for code, shift, phase in zip(coding.codes(chirps), coding.shifts(fft_size), coding.initial_phase):
        measured = np.fft.fft(sequence * code, n=fft_size)
        expected = np.exp(1j * phase) * np.roll(base_spectrum, int(shift))
        shift_error = max(shift_error, float(np.max(np.abs(measured - expected)) / np.max(np.abs(expected))))
    assert shift_error < 1e-10, shift_error

    results = {"modulation_vs_cyclic_shift_relative_error": shift_error}
    speed_period = cfg.wavelength / (2 * cfg.pri)
    for name, velocity in (("stationary", 0.0), ("moving", 18.0), ("aliased_speed", 70.0)):
        target = Target(40, velocity)
        base = simulate(cfg, [target], chirps, None, 0)
        mixture = encode(base, coding, gains, None, 0)
        _, raw_rd, ranges, velocities = process_ddma(cfg, mixture, fft_size)
        _, reference_rd, _, _ = process_ddma(cfg, base, fft_size)
        aligned, score = align_and_score(raw_rd, coding)
        d, r = np.unravel_index(np.argmax(score), score.shape)
        pred = prediction(cfg, target, chirps)
        range_error = abs(ranges[r] - pred["predicted_apparent_range_m"])
        velocity_error = abs(wrap_velocity(velocities[d] - pred["predicted_aliased_velocity_m_s"], speed_period))
        recovered = aligned[:, d, r]
        reference = gains * reference_rd[d, r]
        channel_error = float(np.linalg.norm(recovered - reference) / np.linalg.norm(reference))
        assert range_error <= cfg.range_step, (name, range_error)
        assert velocity_error <= speed_period / fft_size, (name, velocity_error)
        assert channel_error < 0.02, (name, channel_error)
        results[name] = {"range_error_m": float(range_error), "velocity_error_m_s": float(velocity_error),
                         "recovered_complex_channel_relative_error": channel_error}
    return results


def draw(cfg: Config, coding: DDMA, gains: np.ndarray, rd: np.ndarray,
         reference_rd: np.ndarray, aligned: np.ndarray, score: np.ndarray,
         ranges: np.ndarray, velocities: np.ndarray, targets: list[Target],
         predictions: list[dict], detections: list[tuple[int, int]], output: Path) -> None:
    colors = ["#1768ac", "#e28b28", "#29966b", "#a74793"]
    speed_period = cfg.wavelength / (2 * cfg.pri)
    offsets_v = coding.offsets * speed_period
    range_bins = np.arange(len(ranges))
    doppler_bins = np.arange(len(velocities))

    def to_bin(velocity: float) -> float:
        return float((velocity / speed_period * len(velocities)) % len(velocities))

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5))
    m = np.arange(coding.period + 1)
    for q, color in enumerate(colors):
        degrees = ((coding.starts[q] + coding.steps[q] * m) % coding.levels) * 360 / coding.levels
        axes[0].plot(m, degrees, "o-", color=color, label=f"TX{q + 1}: step={coding.steps[q] * 360 / coding.levels:g} deg")
        axes[1].vlines(coding.offsets[q], 0, 1, color=color, lw=3)
        axes[1].plot(coding.offsets[q], 1, "o", color=color)
        axes[1].annotate(f"TX{q + 1}\n{offsets_v[q]:+.2f} m/s", (coding.offsets[q], 1),
                         xytext=(0, 8), textcoords="offset points", ha="center", color=color)
    axes[0].set(xlabel="chirp index", ylabel="transmit phase (deg)", yticks=[0, 90, 180, 270, 360],
                title=f"Assumption: {coding.levels} phase codes = 360 deg; period = {coding.period} chirps")
    axes[0].legend(ncol=2)
    axes[1].set(xlabel="effective baseband offset (cycles / chirp)", ylabel="code locations", xlim=(-0.56, 0.51),
                ylim=(0, 1.5), yticks=[], title=f"Mixer code sign = {coding.sign:+d}; offsets wrap every 1 cycle/chirp")
    for ax in axes:
        ax.grid(alpha=0.2)
    save(fig, output, "01_ddma_phase_codes.png")

    fig, ax = plt.subplots(figsize=(11, 7))
    mesh = ax.pcolormesh(range_bins, doppler_bins, db(rd), shading="auto", cmap="magma", vmin=-65, vmax=0)
    fig.colorbar(mesh, ax=ax, label="relative amplitude (dB)")
    for i, pred in enumerate(predictions, 1):
        for q, offset in enumerate(offsets_v):
            v = to_bin(pred["predicted_aliased_velocity_m_s"] + offset)
            r = pred["predicted_apparent_range_m"] / cfg.range_step
            ax.plot(r, v, "o", mfc="none", mec="cyan", ms=7)
            ax.annotate(f"T{i}/TX{q + 1}", (r, v), xytext=(6, 5), textcoords="offset points", fontsize=8, color="white")
    ax.set(xlabel="Range Bin", ylabel="Doppler Bin",
           title=f"Raw DDMA 2D FFT: natural bin order; {len(targets)} targets, 4 TX")
    save(fig, output, "02_raw_ddma_range_doppler.png")

    chosen = min(range(len(targets)), key=lambda i: abs(targets[i].velocity_m_s))
    pred = predictions[chosen]
    r_bin = int(np.argmin(np.abs(ranges - pred["predicted_apparent_range_m"])))
    reference_amp = float(np.max(np.abs(rd[:, r_bin])))
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    axes[0].plot(doppler_bins, db(rd[:, r_bin], reference_amp), color="#222222")
    for q, offset in enumerate(offsets_v):
        predicted_v = to_bin(pred["predicted_aliased_velocity_m_s"] + offset)
        axes[0].axvline(predicted_v, color=colors[q], linestyle=":", label=f"TX{q + 1} predicted replica")
        axes[1].plot(doppler_bins, db(aligned[q, :, r_bin], reference_amp), alpha=0.6, color=colors[q],
                     label=f"TX{q + 1} code offset removed")
    axes[1].plot(doppler_bins, db(np.sqrt(score[:, r_bin]), reference_amp), color="black", lw=2,
                 label="four-location consistency score")
    axes[1].axvline(to_bin(pred["predicted_aliased_velocity_m_s"]), color="red", linestyle="--", label="expected physical Doppler bin")
    axes[0].set_title(f"Raw Doppler at range bin {r_bin}: one target, four coded peaks")
    axes[1].set_title("Shifting each branch leaves other replicas; the score tests joint alignment")
    axes[1].set_xlabel("Doppler Bin (natural FFT order)")
    for ax in axes:
        ax.set(ylabel="relative response (dB)", ylim=(-75, 5))
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2, loc="lower right")
    save(fig, output, "03_pattern_alignment.png")

    fig, ax = plt.subplots(figsize=(11, 7))
    mesh = ax.pcolormesh(range_bins, doppler_bins, db(np.sqrt(score)), shading="auto", cmap="magma", vmin=-65, vmax=0)
    fig.colorbar(mesh, ax=ax, label="normalized pattern score (dB)")
    for i, (target, pred) in enumerate(zip(targets, predictions), 1):
        r = pred["predicted_apparent_range_m"] / cfg.range_step
        v = to_bin(pred["predicted_aliased_velocity_m_s"])
        ax.plot(r, v, "o", mfc="none", mec="cyan", ms=10)
        ax.annotate(f"T{i}: actual {target.velocity_m_s:g} m/s", (r, v), xytext=(7, 10),
                    textcoords="offset points", fontsize=9, color="white")
    if detections:
        ax.scatter([r for d, r in detections], [d for d, r in detections],
                   marker="x", s=35, color="lime", label="threshold + peak suppression")
        ax.legend(loc="upper right")
    ax.set(xlabel="Range Bin", ylabel="Doppler Bin",
           title="DDMA pattern-consistency detection score: natural bin order")
    save(fig, output, "04_ddma_decoded_map.png")

    if detections:
        d, r = detections[0]
        # 参考信号仅在这个验证图里使用；不进入图样评分或检测。
        recovered_gains = aligned[:, d, r] / reference_rd[d, r]
        q = np.arange(1, 5)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        axes[0].plot(q, np.abs(gains), "o-", label="simulated reference")
        axes[0].plot(q, np.abs(recovered_gains), "x--", ms=9, label="recovered / uncoded reference")
        axes[1].plot(q, np.angle(gains, deg=True), "o-", label="simulated reference")
        axes[1].plot(q, np.angle(recovered_gains, deg=True), "x--", ms=9, label="initial coding phase removed")
        axes[0].set_ylabel("relative channel magnitude")
        axes[1].set_ylabel("channel phase (deg)")
        for ax in axes:
            ax.set(xlabel="TX index (one RX)", xticks=q)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=9)
        fig.suptitle("Validation only: complex TX-channel readout at the strongest detection")
        save(fig, output, "05_virtual_channel_validation.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=PUBLIC_CONFIG)
    parser.add_argument("--output", type=Path, default=None, help="Custom configurations must write outside the repository.")
    parser.add_argument("--chirps", type=int, default=None)
    parser.add_argument("--doppler-fft", type=int, default=None)
    parser.add_argument("--phase-levels", type=int, default=256, help="Assumed code modulus for a full 2pi turn.")
    parser.add_argument("--code-sign", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--snr-db", type=float, default=15)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--threshold-db", type=float, default=-25, help="Relative pattern-power threshold; not CFAR.")
    parser.add_argument("--target", action="append", type=target_arg)
    args = parser.parse_args()
    cfg = Config.read(args.config)
    args.chirps = cfg.doppler_fft if args.chirps is None else args.chirps
    try:
        args.output = choose_output(args.config, args.output, "ddma_results")
    except ValueError as exc:
        parser.error(str(exc))
    header = args.config.read_text(encoding="utf-8-sig")
    coding = DDMA(read_codes(header, "tx_phase_start"), read_codes(header, "tx_phase_step"), args.phase_levels, args.code_sign)
    if len(coding.starts) != 4 or len(coding.steps) != 4:
        parser.error("This demonstration expects four TX phase entries.")
    if coding.levels <= 0 or np.any(coding.starts < 0) or np.any(coding.steps < 0) or np.any(coding.starts >= coding.levels) or np.any(coding.steps >= coding.levels):
        parser.error("Phase codes must be in [0, phase-levels).")
    fft_size = cfg.doppler_fft if args.doppler_fft is None else args.doppler_fft
    if not 64 <= args.chirps <= 4096 or fft_size < args.chirps or args.chirps % coding.period or fft_size % coding.period:
        parser.error(f"Use 64..4096 chirps; actual count and FFT length must be multiples of {coding.period}; FFT >= actual count.")
    if len(set(coding.shifts(fft_size).tolist())) != 4:
        parser.error("Two TX codes have identical Doppler offsets and cannot be separated by this method.")
    symmetry = coding.rotational_symmetries(fft_size)
    if symmetry != [0]:
        parser.error("This full-band pattern decoder needs asymmetric code locations; this pattern has rotational ambiguity.")
    if not args.threshold_db < 0:
        parser.error("Use a negative relative detection threshold in dB.")
    targets = args.target or [Target(20, -12, 1.0, 0.2), Target(45, 18, 0.8, 1.0), Target(70, 70, 0.6, 2.0)]
    if any(not np.isfinite([t.range_m, t.velocity_m_s, t.amplitude]).all() or t.amplitude <= 0 for t in targets):
        parser.error("Targets must have finite values and positive amplitudes.")
    predictions = [prediction(cfg, target, args.chirps) for target in targets]
    if any(not 0 <= p["predicted_apparent_range_m"] < (cfg.range_bins - 1) * cfg.range_step for p in predictions):
        parser.error("Target range peaks must stay within the retained range bins.")
    # 任意已知通道系数，只用于教学和验证；不代表真实天线阵列/角度。
    gains = np.array([1.0, 0.9, 0.8, 0.7]) * np.exp(1j * np.array([0.0, 0.35, -0.55, 0.8]))
    base = simulate(cfg, targets, args.chirps, None, args.seed)
    mixture = encode(base, coding, gains, args.snr_db, args.seed)
    _, raw_rd, ranges, velocities = process_ddma(cfg, mixture, fft_size)
    aligned, score = align_and_score(raw_rd, coding)
    detections = detect(score, args.threshold_db)
    # 到此检测已经完成，下面才计算仅供验证和画图使用的无编码参考。
    _, reference_rd, _, _ = process_ddma(cfg, base, fft_size)
    checks = numerical_checks(cfg, coding, args.chirps, fft_size, gains)
    args.output.mkdir(parents=True, exist_ok=True)
    draw(cfg, coding, gains, raw_rd, reference_rd, aligned, score, ranges, velocities,
         targets, predictions, detections, args.output)

    speed_period = cfg.wavelength / (2 * cfg.pri)
    channels = []
    for q in range(4):
        channels.append({"tx": q + 1, "phase_start_code": int(coding.starts[q]), "phase_step_code": int(coding.steps[q]),
                         "phase_start_deg": float(coding.starts[q] * 360 / coding.levels),
                         "phase_step_deg": float(coding.steps[q] * 360 / coding.levels),
                         "effective_offset_cycles_per_chirp": float(coding.offsets[q]),
                         "doppler_shift_bins": int(coding.shifts(fft_size)[q]),
                         "doppler_shift_bins_modulo_fft": int(coding.shifts(fft_size)[q] % fft_size),
                         "equivalent_velocity_shift_m_s": float(coding.offsets[q] * speed_period)})
    rows = []
    for i, (d, r) in enumerate(detections, 1):
        row = {"detection": i, "range_bin": r, "doppler_bin": d, "range_m": float(ranges[r]),
               "velocity_m_s": float(velocities[d]), "relative_score_db": float(10 * np.log10(score[d, r] / score.max()))}
        rows.append(row)
    with (args.output / "detections.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["detection", "range_bin", "doppler_bin", "range_m", "velocity_m_s", "relative_score_db"])
        writer.writeheader()
        writer.writerows(rows)
    with (args.output / "virtual_channels.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["detection", "tx", "rx", "real", "imag", "magnitude", "phase_deg"])
        for i, (d, r) in enumerate(detections, 1):
            for q, value in enumerate(aligned[:, d, r], 1):
                writer.writerow([i, q, 1, value.real, value.imag, abs(value), np.angle(value, deg=True)])
    report = {
        "model": "4 simultaneous TX, 1 equivalent RX; continuous per-chirp phase ramps",
        "fft_output_order": "natural, unshifted",
        "heatmap_axes": {"x": "Range Bin", "y": "Doppler Bin"},
        "phase_levels_assumed": coding.levels, "code_sign": coding.sign, "code_period_chirps": coding.period,
        "actual_chirps_assumed": args.chirps, "doppler_fft": fft_size, "channels": channels,
        "rotational_symmetries_bins": symmetry, "coherent_time_ms": args.chirps * cfg.pri * 1e3,
        "nominal_velocity_resolution_m_s": speed_period / args.chirps,
        "full_PRF_velocity_half_range_m_s": speed_period / 2,
        "detection_method": "minimum power across four aligned code locations, relative threshold, peak suppression",
        "uses_target_truth_in_detector": False, "detections": rows, "predictions_for_validation_only": predictions,
        "numerical_checks": checks,
        "limitations": ["The bundled configuration and all displayed example parameters are entirely synthetic.",
                        "Four simultaneous TX and phase code units are teaching assumptions.",
                        "Sparse resolvable replicas are required; same-range multi-target collisions can create ambiguity.",
                        "Alignment alone does not isolate the full TX time series; complex channels are read only at detections.",
                        "The code pattern resolves TX-label ambiguity, not aliases beyond the original PRF interval.",
                        "No DOA, actual antenna geometry, RF filters, fixed-point arithmetic or hardware CFAR is emulated."],
        "sources": ["https://www.mathworks.com/help/radar/ug/simulate-an-automotive-4d-imaging-mimo-radar.html"]}
    (args.output / "ddma_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = ["# DDMA 本次结果", "", "4TX 同时发射、1RX 等效接收；仓库自带波形和相位码均为虚构教学值。", "",
             "## 相位码与频移", "", "| TX | 初相位 ° | 步进 ° | 有效 FFT 搬移格点 | 等效速度搬移 m/s |", "|---|---:|---:|---:|---:|"]
    lines += [f"| TX{ch['tx']} | {ch['phase_start_deg']:.4f} | {ch['phase_step_deg']:g} | {ch['doppler_shift_bins']} | {ch['equivalent_velocity_shift_m_s']:.4f} |" for ch in channels]
    lines += ["", "## 检测读数", "", "| 编号 | Range Bin | Doppler Bin | 表观距离 m | PRF 周期内速度 m/s | 相对图样得分 dB |", "|---|---:|---:|---:|---:|---:|"]
    lines += [f"| {row['detection']} | {row['range_bin']} | {row['doppler_bin']} | {row['range_m']:.3f} | {row['velocity_m_s']:.3f} | {row['relative_score_db']:.2f} |" for row in rows]
    lines += ["", "检测使用图样相对阈值，不是硬件 CFAR；读数不依赖真值目标数。完整说明见项目中的 DDMA_GUIDE.md。", ""]
    for filename in ("01_ddma_phase_codes", "02_raw_ddma_range_doppler", "03_pattern_alignment", "04_ddma_decoded_map", "05_virtual_channel_validation"):
        if (args.output / f"{filename}.png").exists():
            lines += [f"![{filename}]({filename}.png)", ""]
    (args.output / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"channels": channels, "detections": rows, "checks": checks}, ensure_ascii=False, indent=2))
    print(f"DDMA results saved to {args.output.resolve()}")


if __name__ == "__main__":
    main()

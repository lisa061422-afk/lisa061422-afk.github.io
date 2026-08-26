"""Reproduce scalar nonlinear timing-perturbation examples for the website.

The numerical integration follows the exploratory MATLAB scripts: explicit
Euler integration with dt = 1e-4 s, piecewise-constant nominal controls, and
delayed control updates.  At each delayed update, the scalar Jacobian-based
LTV gain is recomputed along the nominal trajectory.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from PIL import Image, ImageDraw, ImageFont


DT = 1e-4
TF = 10.0
B = 1.0
UPDATE_TIMES = np.array([1.0, 2.4, 3.8, 5.2, 6.6, 8.0, 9.4, 10.0])


@dataclass(frozen=True)
class Case:
    slug: str
    title: str
    subtitle: str
    model_html: str
    x0: float
    controls: np.ndarray
    delays: np.ndarray
    drift: Callable[[float, float], float]
    jacobian: Callable[[np.ndarray, np.ndarray], np.ndarray]


def control_index(t: np.ndarray | float, update_times: np.ndarray) -> np.ndarray:
    """Match MATLAB's `if t > Time(j)` update convention."""
    return np.searchsorted(update_times, t, side="left")


def simulate_nominal(case: Case, time: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    indices = np.minimum(control_index(time, UPDATE_TIMES), len(case.controls) - 1)
    control = case.controls[indices]
    state = np.empty_like(time)
    state[0] = case.x0
    for i in range(len(time) - 1):
        state[i + 1] = state[i] + DT * (case.drift(time[i], state[i]) + B * control[i])
    return state, control


def scalar_ltv_gain(case: Case, nominal_state: np.ndarray, tau: float, target: float) -> float:
    start = int(round(tau / DT))
    stop = int(round(target / DT))
    if stop <= start:
        raise ValueError(f"Delay leaves no correction interval: tau={tau}, target={target}")

    interval_time = np.arange(start, stop + 1) * DT
    a_values = case.jacobian(interval_time, nominal_state[start : stop + 1])
    cumulative = np.zeros_like(a_values)
    cumulative[1:] = np.cumsum(0.5 * (a_values[:-1] + a_values[1:]) * DT)
    kernel = np.exp(-cumulative)
    integral = np.trapezoid(kernel, dx=DT)
    return 1.0 / (B * integral)


def simulate_delayed(
    case: Case,
    time: np.ndarray,
    nominal_state: np.ndarray,
    corrected: bool,
) -> tuple[np.ndarray, np.ndarray, list[float]]:
    actual_updates = UPDATE_TIMES[:-1] + case.delays
    event_indices = np.rint(actual_updates / DT).astype(int)
    gains = [
        scalar_ltv_gain(case, nominal_state, actual_updates[k], UPDATE_TIMES[k + 1])
        for k in range(len(actual_updates))
    ]

    state = np.empty_like(time)
    control = np.empty_like(time)
    state[0] = case.x0
    current_control = float(case.controls[0])
    next_event = 0

    for i in range(len(time) - 1):
        if next_event < len(event_indices) and i >= event_indices[next_event]:
            current_control = float(case.controls[next_event + 1])
            if corrected:
                error = state[i] - nominal_state[i]
                current_control -= gains[next_event] * error
            next_event += 1
        control[i] = current_control
        state[i + 1] = state[i] + DT * (case.drift(time[i], state[i]) + B * current_control)

    control[-1] = current_control
    return state, control, gains


def nice_ticks(low: float, high: float, count: int = 7) -> np.ndarray:
    span = max(high - low, 1e-9)
    raw_step = span / max(count - 1, 1)
    magnitude = 10 ** math.floor(math.log10(raw_step))
    residual = raw_step / magnitude
    if residual <= 1:
        step = magnitude
    elif residual <= 2:
        step = 2 * magnitude
    elif residual <= 5:
        step = 5 * magnitude
    else:
        step = 10 * magnitude
    start = math.floor(low / step) * step
    end = math.ceil(high / step) * step
    return np.arange(start, end + 0.5 * step, step)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    filename = "arialbd.ttf" if bold else "arial.ttf"
    path = Path("C:/Windows/Fonts") / filename
    return ImageFont.truetype(str(path), size=size)


def draw_patterned_polyline(
    draw: ImageDraw.ImageDraw,
    points: list[tuple[float, float]],
    fill: tuple[int, int, int],
    width: int,
    pattern: tuple[float, ...],
) -> None:
    pattern_index = 0
    remaining = pattern[0]
    drawing = True
    for p0, p1 in zip(points[:-1], points[1:]):
        x0, y0 = p0
        x1, y1 = p1
        dx, dy = x1 - x0, y1 - y0
        segment_length = math.hypot(dx, dy)
        if segment_length == 0:
            continue
        consumed = 0.0
        while consumed < segment_length:
            take = min(remaining, segment_length - consumed)
            a = consumed / segment_length
            b = (consumed + take) / segment_length
            start = (x0 + a * dx, y0 + a * dy)
            end = (x0 + b * dx, y0 + b * dy)
            if drawing:
                draw.line([start, end], fill=fill, width=width)
            consumed += take
            remaining -= take
            if remaining <= 1e-8:
                pattern_index = (pattern_index + 1) % len(pattern)
                remaining = pattern[pattern_index]
                drawing = pattern_index % 2 == 0


def format_tick(value: float) -> str:
    if abs(value) < 1e-10:
        return "0"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.1f}"


def render_case(
    case: Case,
    time: np.ndarray,
    nominal_state: np.ndarray,
    nominal_control: np.ndarray,
    corrected_state: np.ndarray,
    corrected_control: np.ndarray,
    output_path: Path,
) -> None:
    width, height = 1800, 1160
    plot_left, plot_right = 175, width - 70
    state_top, state_bottom = 145, 700
    control_top, control_bottom = 790, 990
    sampled = slice(None, None, 10)

    state_values = np.concatenate([nominal_state[sampled], corrected_state[sampled]])
    control_values = np.concatenate([nominal_control[sampled], corrected_control[sampled]])
    state_pad = max(0.08 * np.ptp(state_values), 0.12)
    control_pad = max(0.10 * np.ptp(control_values), 0.08)
    state_ticks = nice_ticks(float(np.min(state_values) - state_pad), float(np.max(state_values) + state_pad))
    control_ticks = nice_ticks(float(np.min(control_values) - control_pad), float(np.max(control_values) + control_pad), count=5)
    state_low, state_high = float(state_ticks[0]), float(state_ticks[-1])
    control_low, control_high = float(control_ticks[0]), float(control_ticks[-1])

    def px_x(values: np.ndarray | float) -> np.ndarray | float:
        return plot_left + np.asarray(values) / TF * (plot_right - plot_left)

    def state_y(values: np.ndarray | float) -> np.ndarray | float:
        return state_bottom - (np.asarray(values) - state_low) / (state_high - state_low) * (state_bottom - state_top)

    def control_y(values: np.ndarray | float) -> np.ndarray | float:
        return control_bottom - (np.asarray(values) - control_low) / (control_high - control_low) * (control_bottom - control_top)

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    grid = (216, 219, 223)
    axis = (32, 37, 42)
    text = (34, 39, 44)
    muted = (91, 99, 110)
    blue = (74, 78, 255)
    red = (245, 20, 18)
    black = (12, 14, 16)

    x_ticks = np.arange(0.0, TF + 0.01, 2.0)
    for tick in x_ticks:
        x = float(px_x(tick))
        draw.line([(x, state_top), (x, state_bottom)], fill=grid, width=2)
        draw.line([(x, control_top), (x, control_bottom)], fill=grid, width=2)
        label = format_tick(float(tick))
        bbox = draw.textbbox((0, 0), label, font=font(26))
        draw.text((x - (bbox[2] - bbox[0]) / 2, control_bottom + 16), label, fill=text, font=font(26))

    for tick in state_ticks:
        y = float(state_y(tick))
        draw.line([(plot_left, y), (plot_right, y)], fill=grid, width=2)
        label = format_tick(float(tick))
        bbox = draw.textbbox((0, 0), label, font=font(25))
        draw.text((plot_left - 18 - (bbox[2] - bbox[0]), y - 15), label, fill=text, font=font(25))

    for tick in control_ticks:
        y = float(control_y(tick))
        draw.line([(plot_left, y), (plot_right, y)], fill=grid, width=2)
        label = format_tick(float(tick))
        bbox = draw.textbbox((0, 0), label, font=font(23))
        draw.text((plot_left - 18 - (bbox[2] - bbox[0]), y - 14), label, fill=text, font=font(23))

    draw.rectangle([(plot_left, state_top), (plot_right, state_bottom)], outline=axis, width=3)
    draw.rectangle([(plot_left, control_top), (plot_right, control_bottom)], outline=axis, width=3)

    for update_time in UPDATE_TIMES:
        x = float(px_x(update_time))
        draw_patterned_polyline(draw, [(x, state_top), (x, state_bottom)], blue, 4, (25, 15))
        draw_patterned_polyline(draw, [(x, control_top), (x, control_bottom)], blue, 4, (25, 15))

    t_sample = time[sampled]
    nominal_points = list(zip(px_x(t_sample).tolist(), state_y(nominal_state[sampled]).tolist()))
    corrected_points = list(zip(px_x(t_sample).tolist(), state_y(corrected_state[sampled]).tolist()))
    nominal_u_points = list(zip(px_x(t_sample).tolist(), control_y(nominal_control[sampled]).tolist()))
    corrected_u_points = list(zip(px_x(t_sample).tolist(), control_y(corrected_control[sampled]).tolist()))
    draw.line(nominal_points, fill=black, width=6, joint="curve")
    draw.line(corrected_points, fill=red, width=6, joint="curve")
    draw_patterned_polyline(draw, nominal_u_points, black, 5, (24, 9, 4, 9))
    draw_patterned_polyline(draw, corrected_u_points, red, 5, (24, 9, 4, 9))

    title_box = draw.textbbox((0, 0), case.title, font=font(38, bold=True))
    draw.text(((width - (title_box[2] - title_box[0])) / 2, 26), case.title, fill=text, font=font(38, bold=True))
    subtitle_box = draw.textbbox((0, 0), case.subtitle, font=font(24))
    draw.text(((width - (subtitle_box[2] - subtitle_box[0])) / 2, 78), case.subtitle, fill=muted, font=font(24))

    x_label = "Time (s)"
    x_box = draw.textbbox((0, 0), x_label, font=font(34))
    draw.text(((width - (x_box[2] - x_box[0])) / 2, height - 75), x_label, fill=(32, 96, 46), font=font(34))

    for label, center_y in [("State", (state_top + state_bottom) // 2), ("Control", (control_top + control_bottom) // 2)]:
        layer = Image.new("RGBA", (220, 65), (255, 255, 255, 0))
        layer_draw = ImageDraw.Draw(layer)
        layer_draw.text((0, 4), label, fill=text, font=font(31))
        layer = layer.rotate(90, expand=True)
        image.paste(layer, (35, int(center_y - layer.height / 2)), layer)
    draw = ImageDraw.Draw(image)

    state_legend = [
        ("Nominal state", black, "solid"),
        ("Corrected delayed state", red, "solid"),
        ("Nominal update time", blue, "dashed"),
    ]
    legend_width, legend_height = 435, 132
    legend_x, legend_y = plot_right - legend_width - 20, state_top + 18
    draw.rounded_rectangle(
        [(legend_x, legend_y), (legend_x + legend_width, legend_y + legend_height)],
        radius=9, fill=(255, 255, 255), outline=(154, 160, 168), width=2,
    )
    for row, (label, color, style) in enumerate(state_legend):
        y = legend_y + 25 + row * 38
        p0, p1 = (legend_x + 20, y), (legend_x + 105, y)
        if style == "solid":
            draw.line([p0, p1], fill=color, width=5)
        else:
            draw_patterned_polyline(draw, [p0, p1], color, 5, (20, 12))
        draw.text((legend_x + 122, y - 14), label, fill=text, font=font(22))

    control_legend = [
        ("Nominal control", black),
        ("Corrected delayed control", red),
    ]
    control_legend_width, control_legend_height = 435, 92
    control_legend_x = plot_right - control_legend_width - 20
    control_legend_y = control_top + 18
    draw.rounded_rectangle(
        [(control_legend_x, control_legend_y), (control_legend_x + control_legend_width, control_legend_y + control_legend_height)],
        radius=9, fill=(255, 255, 255), outline=(154, 160, 168), width=2,
    )
    for row, (label, color) in enumerate(control_legend):
        y = control_legend_y + 25 + row * 38
        p0, p1 = (control_legend_x + 20, y), (control_legend_x + 105, y)
        draw_patterned_polyline(draw, [p0, p1], color, 5, (20, 8, 4, 8))
        draw.text((control_legend_x + 122, y - 14), label, fill=text, font=font(22))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(output_path, format="PNG", optimize=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    cases = [
        Case(
            slug="forced-pendulum",
            title="Case 1: Periodically Forced Nonlinear Drift",
            subtitle="xdot = -0.70 sin(x) + 0.55 sin(1.25t) + 0.18 sin(3.2t) + u",
            model_html=r"\(\dot x=-0.70\sin x+0.55\sin(1.25t)+0.18\sin(3.2t)+u\), \(A(t)=-0.70\cos(x^n(t))\)",
            x0=0.90,
            controls=np.array([0.0, 0.18, -0.12, 0.20, -0.15, 0.12, -0.10, 0.0]),
            delays=np.array([0.18, 0.00, 0.22, 0.00, 0.16, 0.00, 0.20]),
            drift=lambda t, x: -0.70 * math.sin(x) + 0.55 * math.sin(1.25 * t) + 0.18 * math.sin(3.2 * t),
            jacobian=lambda t, x: -0.70 * np.cos(x),
        ),
        Case(
            slug="forced-bistable-cubic",
            title="Case 2: Periodically Forced Bistable Cubic System",
            subtitle="xdot = 0.75x - 0.35x^3 + 0.45 cos(1.1t) + 0.16 sin(2.7t) + u",
            model_html=r"\(\dot x=0.75x-0.35x^3+0.45\cos(1.1t)+0.16\sin(2.7t)+u\), \(A(t)=0.75-1.05(x^n(t))^2\)",
            x0=-0.80,
            controls=np.array([0.0, -0.10, 0.12, -0.08, 0.10, -0.12, 0.08, 0.0]),
            delays=np.array([0.15, 0.00, 0.20, 0.00, 0.18, 0.00, 0.16]),
            drift=lambda t, x: 0.75 * x - 0.35 * x**3 + 0.45 * math.cos(1.1 * t) + 0.16 * math.sin(2.7 * t),
            jacobian=lambda t, x: 0.75 - 1.05 * x**2,
        ),
        Case(
            slug="modulated-multiharmonic",
            title="Case 3: Time-Modulated Multi-Harmonic Drift",
            subtitle="xdot = [0.45 + 0.28 sin(0.85t)] sin(2x) + periodic forcing + u",
            model_html=r"\(\dot x=[0.45+0.28\sin(0.85t)]\sin(2x)+0.32\cos(1.75t)+0.14\sin(3.6t)+u\)",
            x0=0.20,
            controls=np.array([0.0, 0.12, -0.10, 0.14, -0.12, 0.10, -0.08, 0.0]),
            delays=np.array([0.20, 0.00, 0.18, 0.00, 0.22, 0.00, 0.16]),
            drift=lambda t, x: (0.45 + 0.28 * math.sin(0.85 * t)) * math.sin(2.0 * x) + 0.32 * math.cos(1.75 * t) + 0.14 * math.sin(3.6 * t),
            jacobian=lambda t, x: 2.0 * (0.45 + 0.28 * np.sin(0.85 * t)) * np.cos(2.0 * x),
        ),
    ]

    time = np.arange(0.0, TF + 0.5 * DT, DT)
    summaries = []

    for case in cases:
        nominal_state, nominal_control = simulate_nominal(case, time)
        corrected_state, corrected_control, gains = simulate_delayed(
            case, time, nominal_state, corrected=True
        )
        uncorrected_state, _, _ = simulate_delayed(case, time, nominal_state, corrected=False)

        corrected_error = corrected_state - nominal_state
        uncorrected_error = uncorrected_state - nominal_state
        corrected_rms = float(np.sqrt(np.mean(corrected_error**2)))
        uncorrected_rms = float(np.sqrt(np.mean(uncorrected_error**2)))
        improvement = 100.0 * (1.0 - corrected_rms / uncorrected_rms)

        output_path = args.output_dir / f"sim_{case.slug}.png"
        render_case(
            case,
            time,
            nominal_state,
            nominal_control,
            corrected_state,
            corrected_control,
            output_path,
        )

        summary = {
            "slug": case.slug,
            "title": case.title,
            "model_html": case.model_html,
            "delays_seconds": case.delays.tolist(),
            "gains": [round(value, 6) for value in gains],
            "rms_error_uncorrected": uncorrected_rms,
            "rms_error_corrected": corrected_rms,
            "rms_error_reduction_percent": improvement,
            "max_abs_error_corrected": float(np.max(np.abs(corrected_error))),
            "final_abs_error_corrected": float(abs(corrected_error[-1])),
            "figure": output_path.name,
        }
        summaries.append(summary)
        print(
            f"{case.slug}: RMS {uncorrected_rms:.5f} -> {corrected_rms:.5f} "
            f"({improvement:.1f}% reduction), max corrected error "
            f"{summary['max_abs_error_corrected']:.5f}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "simulation_summary.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()

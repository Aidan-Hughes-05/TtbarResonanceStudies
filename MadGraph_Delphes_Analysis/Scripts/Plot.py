#!/usr/bin/env python3
"""Plotting script using uproot and matplotlib based on YAML configuration."""

import argparse
import logging
import os
import warnings
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import uproot
import yaml

# Configure logger
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("PlottingScript")


def load_config(config_path: str | Path) -> Dict[str, Any]:
    """Load configuration from a YAML file."""
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def load_tree_data(
    file_path: str | Path,
    tree_name: str,
    observable: str,
    weight_branch: str,
    params: Optional[Dict[str, Any]] = None,
    apply_param_cuts: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Read observable and weight arrays from a ROOT TTree applying parameter filters."""
    branches_to_read = [observable, weight_branch]
    filter_params = params or {}

    with uproot.open(file_path) as root_file:
        if tree_name not in root_file:
            raise KeyError(f"TTree '{tree_name}' not found in {file_path}")
        tree = root_file[tree_name]
        available_keys = set(tree.keys())

        # Check required observable and weight branches
        for b in branches_to_read:
            if b not in available_keys:
                raise KeyError(f"Branch '{b}' not found in TTree '{tree_name}' of {file_path}")

        # Include filter parameters present in this tree if cuts enabled
        active_cuts = {}
        if apply_param_cuts:
            for param_key, param_val in filter_params.items():
                if param_key in available_keys:
                    active_cuts[param_key] = param_val
                    if param_key not in branches_to_read:
                        branches_to_read.append(param_key)

        arrays = tree.arrays(branches_to_read, library="np")

    # Apply cuts
    mask = np.ones(len(arrays[observable]), dtype=bool)
    if apply_param_cuts and active_cuts:
        for k, v in active_cuts.items():
            val_arr = arrays[k]
            if np.issubdtype(val_arr.dtype, np.floating):
                mask = mask & np.isclose(val_arr, float(v))
            else:
                mask = mask & (val_arr == v)

    values = arrays[observable][mask]
    weights = arrays[weight_branch][mask]

    return values, weights


def compute_histogram(
    values: np.ndarray,
    weights: np.ndarray,
    n_bins: int,
    x_range: Tuple[float, float],
    is_density: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute binned counts and optionally normalize by the absolute sum of bin contents * bin width."""
    counts, bin_edges = np.histogram(
        values, bins=n_bins, range=x_range, weights=weights
    )

    if is_density:
        bin_widths = np.diff(bin_edges)
        # Area defined as abs(negative + positive) contributions
        abs_area = np.sum(np.abs(counts) * bin_widths)
        if abs_area > 0:
            counts = counts / abs_area
        else:
            logger.warning("Total area is 0; skipping density normalization.")

    return counts, bin_edges


def format_param_title(params: Optional[Dict[str, Any]]) -> str:
    """Format the first two parameter values using the configured display labels."""
    if not params:
        return ""

    parameter_items = [(key, value) for key, value in params.items() if key != "labels"][:2]
    labels = params.get("labels", [])
    parts = []
    for index, (key, value) in enumerate(parameter_items):
        display_label = labels[index] if index < len(labels) else key
        key_str = f"${display_label}$" if any(char in str(display_label) for char in "_^") else str(display_label)
        parts.append(f"{key_str} = {value}")
    return ", ".join(parts)


def format_param_filename(params: Optional[Dict[str, Any]]) -> str:
    """Format parameter dictionary into a file-system safe string."""
    if not params:
        return "inclusive"
    parts = []
    for k, v in sorted(params.items()):
        if k == "labels":
            continue
        val_str = str(v).replace(".", "p")
        parts.append(f"{k}_{val_str}")
    return "_".join(parts)


def plot_observable(
    config: Dict[str, Any],
    obs_config: Dict[str, Any],
    param_dict: Optional[Dict[str, Any]],
    output_dir: Path,
) -> None:
    """Generate and save the scientific plot for a given observable and parameter point."""
    obs_name = obs_config["name"]
    x_range = tuple(param_dict.get("x_range", obs_config.get("x_range", None)))
    if x_range[0] is None or x_range[1] is None:
        raise ValueError(f"x_range must be specified for observable '{obs_name}' or in parameters.")
    x_label = obs_config.get("x_label", obs_name)
    y_label = obs_config.get("y_label", "Events")
    weights_branch = obs_config.get("weights", "lumiWeight")
    n_bins = obs_config.get("n_bins", 50)
    font_size = obs_config.get("font_size", 16)
    is_density = (y_label.strip().lower() == "density")

    tree_name = config.get("ttree", "Ttbar")
    bg_path = config["background_path"]
    sig_path = config["signal_path"]
    int_path = config["interference_path"]
    skip_background = bool(param_dict.get("skip_background", False))

    # 1. Load data for components
    # Background does not depend on resonance parameters (apply_param_cuts=False)
    if skip_background:
        bg_vals = np.array([])
        bg_weights = np.array([])
    else:
        bg_vals, bg_weights = load_tree_data(
            bg_path, tree_name, obs_name, weights_branch, param_dict, apply_param_cuts=False
        )
    sig_vals, sig_weights = load_tree_data(
        sig_path, tree_name, obs_name, weights_branch, param_dict, apply_param_cuts=True
    )
    int_vals, int_weights = load_tree_data(
        int_path, tree_name, obs_name, weights_branch, param_dict, apply_param_cuts=True
    )

    # 2. Combine Signal + Interference
    sig_int_vals = np.concatenate([sig_vals, int_vals])
    sig_int_weights = np.concatenate([sig_weights, int_weights])

    # 3. Setup styling
    plt.rcParams.update({
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size - 2,
        "ytick.labelsize": font_size - 2,
        "legend.fontsize": font_size - 4,
    })

    fig, ax = plt.subplots(figsize=(10, 7), dpi=150)

    # Color palette
    colors = {
        "Background": "#1f77b4",  # Blue
        "Signal": "#2ca02c",      # Green
        "Interference": "#d62728", # Red
        "Signal + Interference": "#9467bd",  # Purple
    }

    # Line styles
    ls_pos = "--"
    ls_neg = ":"
    lw_total = 2.0
    lw_sub = 1.4

    # Data components dictionary
    components = {
        "Signal": (sig_vals, sig_weights),
        "Interference": (int_vals, int_weights),
    }
    if not skip_background:
        components["Background"] = (bg_vals, bg_weights)

    # Plot total components
    for comp_name, (vals, w) in components.items():
        if len(vals) == 0:
            logger.warning(f"No events found for {comp_name} at {param_dict}")
            continue
        counts, bin_edges = compute_histogram(vals, w, n_bins, x_range, is_density=is_density)
        ax.stairs(
            counts,
            bin_edges,
            label=f"{comp_name}",
            color=colors[comp_name],
            linestyle="-",
            linewidth=lw_total,
        )

    # Plot Signal + Interference total
    if len(sig_int_vals) > 0:
        counts_si, bin_edges_si = compute_histogram(
            sig_int_vals, sig_int_weights, n_bins, x_range, is_density=is_density
        )
        ax.stairs(
            counts_si,
            bin_edges_si,
            label="Signal + Interference",
            color=colors["Signal + Interference"],
            linestyle="-",
            linewidth=lw_total,
        )
    else:
        logger.warning(f"No events found for Signal + Interference at {param_dict}")

    # Plot positive and negative weight components for Background, Signal, Interference
    # Only decompose if a histogram contains BOTH positive and negative weight contributions.
    # If it has only one sign, the solid total histogram already represents it, avoiding redundant duplicate lines.
    for comp_name, (vals, w) in components.items():
        if len(vals) == 0:
            continue

        has_pos = np.any(w > 0)
        has_neg = np.any(w < 0)

        # Decompose only when both signs are present
        if has_pos and has_neg:
            pos_mask = w > 0
            neg_mask = w < 0

            c_pos, be_pos = compute_histogram(
                vals[pos_mask], w[pos_mask], n_bins, x_range, is_density=is_density
            )
            ax.stairs(
                c_pos,
                be_pos,
                label=f"{comp_name} ($w > 0$)",
                color=colors[comp_name],
                linestyle=ls_pos,
                linewidth=lw_sub,
            )

            c_neg, be_neg = compute_histogram(
                vals[neg_mask], w[neg_mask], n_bins, x_range, is_density=is_density
            )
            ax.stairs(
                c_neg,
                be_neg,
                label=f"{comp_name} ($w < 0$)",
                color=colors[comp_name],
                linestyle=ls_neg,
                linewidth=lw_sub,
            )
        else:
            missing_sign = "negative" if not has_neg else "positive"
            logger.warning(
                f"{comp_name} only contains one weight sign (missing {missing_sign} weights) for {param_dict}; "
                f"omitting separate sign-split curves."
            )

    # Configure axes and labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(x_range)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle="-", alpha=0.5)
    ax.grid(True, linestyle="--", alpha=0.4)

    param_title = format_param_title(param_dict)
    if param_title:
        ax.set_title(param_title, pad=12)

    # Place legend outside the main plotting axes so it doesn't cover data
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
        framealpha=0.9,
    )

    fig.tight_layout()

    # Save outputs (both PNG and PDF)
    param_tag = format_param_filename(param_dict)
    out_filename_png = f"{obs_name}_{param_tag}.png"
    out_filename_pdf = f"{obs_name}_{param_tag}.pdf"

    output_path_png = output_dir / out_filename_png
    output_path_pdf = output_dir / out_filename_pdf

    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"Saved plot: {output_path_png} and {output_path_pdf}")


def main() -> None:
    """Main execution entry point."""
    default_config_path = Path(__file__).resolve().parent / "PlotsConfig.yaml"

    parser = argparse.ArgumentParser(description="Create resonance analysis comparison plots.")
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default=str(default_config_path),
        help="Path to YAML configuration file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default=None,
        help="Directory to save generated plots (default: ./plots relative to config directory).",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    config = load_config(config_path)

    if args.output_dir:
        output_dir = Path(args.output_dir).resolve()
    else:
        output_dir = config_path.parent / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Loaded config from {config_path}")
    logger.info(f"Saving output plots to {output_dir}")

    observables = config.get("observables", [])
    if not observables:
        logger.warning("No observables defined in configuration.")
        return

    for obs in observables:
        params_list = obs.get("params", [None])
        if not params_list:
            params_list = [None]

        for param_dict in params_list:
            logger.info(f"Processing observable '{obs.get('name')}' with parameters: {param_dict}")
            plot_observable(config, obs, param_dict, output_dir)


if __name__ == "__main__":
    main()

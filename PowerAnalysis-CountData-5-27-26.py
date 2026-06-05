#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Power analysis for coral embryo count data.

Experimental unit: tube.

Each tube contains a random number of coral between
min_coral_per_tube and max_coral_per_tube.

Each tube contributes one response value:
    proportion viable = viable / total

Power test:
    Welch two-sample t-test comparing tube-level proportions
    between control and treatment groups.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind


def logit(p):
    """Logit transform with clipping for numerical stability."""
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return np.log(p / (1 - p))


def expit(x):
    """Inverse logit."""
    return 1 / (1 + np.exp(-x))


def get_treatment_p(control_p, effect_size, effect_type="absolute"):
    """
    Compute treatment probability from control probability.
    """
    if effect_type == "absolute":
        treatment_p = control_p - effect_size
    elif effect_type == "relative":
        treatment_p = control_p * (1 - effect_size)
    else:
        raise ValueError("effect_type must be 'absolute' or 'relative'")

    return np.clip(treatment_p, 1e-9, 1 - 1e-9)


def simulate_tube_counts(
    n_tubes_per_group,
    min_coral_per_tube,
    max_coral_per_tube,
    control_p=0.80,
    effect_size=0.10,
    effect_type="absolute",
    tube_sd=0.0,
    seed=None
):
    """
    Simulate one experiment.

    One row = one tube.
    Each tube has a random number of coral.
    """
    rng = np.random.default_rng(seed)

    treatment_p = get_treatment_p(
        control_p=control_p,
        effect_size=effect_size,
        effect_type=effect_type
    )

    def make_group(base_p, treatment_label):
        rows = []
        base_logit = logit(base_p)

        for tube_id in range(n_tubes_per_group):

            coral_in_tube = rng.integers(
                min_coral_per_tube,
                max_coral_per_tube + 1
            )

            tube_effect = rng.normal(0, tube_sd)
            p_tube = expit(base_logit + tube_effect)

            viable = rng.binomial(coral_in_tube, p_tube)

            rows.append({
                "treatment": treatment_label,
                "tube_id": tube_id,
                "viable": viable,
                "total": coral_in_tube,
                "prop": viable / coral_in_tube,
                "p_tube": p_tube
            })

        return rows

    control_rows = make_group(control_p, 0)
    treatment_rows = make_group(treatment_p, 1)

    return pd.DataFrame(control_rows + treatment_rows)


def tube_level_test_pvalue(df):
    """
    Compare treatment groups using tube-level proportions.

    Tube is the experimental unit.
    """
    control_props = df[df["treatment"] == 0]["prop"].to_numpy()
    treatment_props = df[df["treatment"] == 1]["prop"].to_numpy()

    if len(control_props) < 2 or len(treatment_props) < 2:
        return 1.0

    if (
        np.var(control_props) == 0
        and np.var(treatment_props) == 0
        and np.mean(control_props) == np.mean(treatment_props)
    ):
        return 1.0

    try:
        _, pval = ttest_ind(
            control_props,
            treatment_props,
            equal_var=False,
            alternative="two-sided"
        )
    except Exception:
        pval = 1.0

    if np.isnan(pval):
        pval = 1.0

    return pval


def estimate_power_counts(
    n_tubes_per_group,
    min_coral_per_tube,
    max_coral_per_tube,
    control_p=0.80,
    effect_size=0.10,
    effect_type="absolute",
    tube_sd=0.0,
    alpha=0.05,
    n_sim=1000,
    seed=123
):
    """
    Estimate statistical power by repeated simulation.
    """
    rng = np.random.default_rng(seed)
    reject = 0

    for _ in range(n_sim):
        sim_seed = rng.integers(0, 1_000_000_000)

        df = simulate_tube_counts(
            n_tubes_per_group=n_tubes_per_group,
            min_coral_per_tube=min_coral_per_tube,
            max_coral_per_tube=max_coral_per_tube,
            control_p=control_p,
            effect_size=effect_size,
            effect_type=effect_type,
            tube_sd=tube_sd,
            seed=sim_seed
        )

        pvalue = tube_level_test_pvalue(df)

        if pvalue < alpha:
            reject += 1

    return reject / n_sim


def find_required_tubes_counts(
    effect_sizes,
    min_coral_per_tube,
    max_coral_per_tube,
    control_p=0.80,
    effect_type="absolute",
    tube_sd=0.0,
    alpha=0.05,
    target_power=0.80,
    n_sim=1000,
    min_tubes=2,
    max_tubes=50,
    seed=123,
    verbose=True
):
    """
    Find the minimum number of tubes per group needed for each effect size.
    """
    results = []

    for effect in effect_sizes:
        found = False

        if verbose:
            if effect_type == "absolute":
                print(f"\nTesting absolute drop of {effect:.2f} ...")
            else:
                print(f"\nTesting relative drop of {effect:.2f} ...")

        for n_tubes in range(min_tubes, max_tubes + 1):
            power = estimate_power_counts(
                n_tubes_per_group=n_tubes,
                min_coral_per_tube=min_coral_per_tube,
                max_coral_per_tube=max_coral_per_tube,
                control_p=control_p,
                effect_size=effect,
                effect_type=effect_type,
                tube_sd=tube_sd,
                alpha=alpha,
                n_sim=n_sim,
                seed=seed + n_tubes + int(effect * 1000)
            )

            if verbose:
                print(f"  tubes/group = {n_tubes:2d}, power = {power:.3f}")

            if power >= target_power:
                treatment_p = get_treatment_p(
                    control_p=control_p,
                    effect_size=effect,
                    effect_type=effect_type
                )

                results.append({
                    "effect_size": effect,
                    "effect_type": effect_type,
                    "control_p": control_p,
                    "treatment_p": treatment_p,
                    "tubes_per_group_needed": n_tubes,
                    "estimated_power": power
                })

                found = True
                break

        if not found:
            treatment_p = get_treatment_p(
                control_p=control_p,
                effect_size=effect,
                effect_type=effect_type
            )

            results.append({
                "effect_size": effect,
                "effect_type": effect_type,
                "control_p": control_p,
                "treatment_p": treatment_p,
                "tubes_per_group_needed": np.nan,
                "estimated_power": np.nan
            })

    return pd.DataFrame(results)


def make_power_curve_counts(
    effect_size,
    min_coral_per_tube,
    max_coral_per_tube,
    control_p=0.80,
    effect_type="absolute",
    tube_sd=0.0,
    alpha=0.05,
    n_sim=1000,
    min_tubes=2,
    max_tubes=50,
    seed=123
):
    """
    Generate a power curve for one effect size.
    """
    curve = []

    for n_tubes in range(min_tubes, max_tubes + 1):
        power = estimate_power_counts(
            n_tubes_per_group=n_tubes,
            min_coral_per_tube=min_coral_per_tube,
            max_coral_per_tube=max_coral_per_tube,
            control_p=control_p,
            effect_size=effect_size,
            effect_type=effect_type,
            tube_sd=tube_sd,
            alpha=alpha,
            n_sim=n_sim,
            seed=seed + n_tubes
        )

        curve.append((n_tubes, power))

    return pd.DataFrame(
        curve,
        columns=["tubes_per_group", "power"]
    )


if __name__ == "__main__":

    # ------------------------------------------------
    # USER SETTINGS
    # ------------------------------------------------

    min_coral_per_tube = 10
    max_coral_per_tube = 10

    control_p = 0.90

    effect_type = "absolute"
    effect_sizes = (0.10, 0.20, 0.30)

    tube_sd = 0.5

    alpha = 0.05
    target_power = 0.80
    n_sim = 1000

    min_tubes = 2
    max_tubes = 50

    # ------------------------------------------------
    # FIND REQUIRED TUBES
    # ------------------------------------------------

    results = find_required_tubes_counts(
        effect_sizes=effect_sizes,
        min_coral_per_tube=min_coral_per_tube,
        max_coral_per_tube=max_coral_per_tube,
        control_p=control_p,
        effect_type=effect_type,
        tube_sd=tube_sd,
        alpha=alpha,
        target_power=target_power,
        n_sim=n_sim,
        min_tubes=min_tubes,
        max_tubes=max_tubes,
        seed=123,
        verbose=True
    )

    print("\nMinimum tubes per group needed:")
    print(results.to_string(index=False))

    # ------------------------------------------------
    # PLOT POWER CURVES
    # ------------------------------------------------

    plt.figure(figsize=(8, 6))

    for effect in effect_sizes:
        curve = make_power_curve_counts(
            effect_size=effect,
            min_coral_per_tube=min_coral_per_tube,
            max_coral_per_tube=max_coral_per_tube,
            control_p=control_p,
            effect_type=effect_type,
            tube_sd=tube_sd,
            alpha=alpha,
            n_sim=n_sim,
            min_tubes=min_tubes,
            max_tubes=max_tubes,
            seed=123
        )

        if effect_type == "absolute":
            label = f"{int(effect * 100)} percentage-point drop"
        else:
            label = f"{int(effect * 100)}% relative drop"

        plt.plot(
            curve["tubes_per_group"],
            curve["power"],
            label=label
        )

    plt.axhline(
        target_power,
        linestyle="--",
        label=f"{int(target_power * 100)}% power"
    )

    plt.xlabel("Tubes per group")
    plt.ylabel("Estimated power")

    plt.title(
        "Tube-level power analysis for coral count data\n"
        f"{min_coral_per_tube}-{max_coral_per_tube} coral/tube, "
        f"control p={control_p:.2f}, "
        f"effect={effect_type}, tube_sd={tube_sd}"
    )

    plt.legend()
    plt.tight_layout()
    plt.show()
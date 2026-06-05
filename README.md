# Coral Embryo Viability Power Analysis

## Overview

This repository contains Python code for conducting simulation-based statistical power analyses for coral embryo viability experiments. The goal is to estimate the number of replicate tubes required per treatment group to achieve a target statistical power (typically 80%) for detecting biologically meaningful reductions in embryo viability.

The analysis treats each tube as the experimental unit and incorporates both binomial sampling variation and among-tube biological variability.

The primary script is:

`PowerAnalysis-CountData-5-27-26.py`

---

## Scientific Background

Coral embryo viability experiments commonly assess the proportion of viable embryos within replicate tubes. Because viability can vary among tubes due to biological and experimental factors, power calculations based solely on individual embryos may substantially overestimate statistical power.

This simulation framework models viability at the tube level, allowing realistic estimation of the number of replicate tubes required to detect specified treatment effects.

---

## Features

* Simulation-based power analysis using Monte Carlo methods.
* Models embryo viability as a binomial process.
* Allows fixed or variable numbers of embryos per tube.
* Incorporates tube-to-tube heterogeneity using random effects on the logit scale.
* Supports both absolute and relative reductions in viability.
* Uses Welch's two-sample t-test on tube-level viability proportions.
* Estimates statistical power across a range of replication levels.
* Generates power curves showing the relationship between replication and statistical power.
* Determines the minimum number of replicate tubes needed to achieve a target power threshold.

---

## Statistical Approach

### Experimental Unit

The experimental unit is the tube.

Each tube contains a user-defined number of embryos bounded by:

* `min_coral_per_tube`
* `max_coral_per_tube`

Each tube contributes a single response value:

```
proportion viable = viable embryos / total embryos
```

### Viability Simulation

For each tube:

1. A baseline control viability probability is specified.
2. A treatment viability probability is generated using either:

   * an absolute reduction, or
   * a relative reduction.
3. Tube-level variation is introduced by adding a normally distributed random effect on the logit scale.
4. Viable embryos are sampled from a binomial distribution.

### Statistical Testing

Treatment and control groups are compared using Welch's two-sample t-test applied to tube-level viability proportions.

Welch's t-test was selected because it:

* treats the tube as the experimental unit,
* does not assume equal variances between groups,
* is robust to unequal variability among treatments.

### Power Estimation

Statistical power is estimated through repeated Monte Carlo simulation.

For each simulated experiment:

1. Control and treatment tubes are generated.
2. A Welch's t-test is performed.
3. The treatment effect is considered detected when:

```
p < 0.05
```

Power is calculated as:

```
Power = Number of significant simulations / Total simulations
```

The procedure is repeated across increasing numbers of replicate tubes to identify the minimum replication needed to achieve the target power.

---

## Requirements

* Python ≥ 3.9
* NumPy
* Pandas
* Matplotlib
* SciPy

Install dependencies:

```bash
pip install numpy pandas matplotlib scipy
```

---

## Usage

Run the script:

```bash
python PowerAnalysis-CountData-5-27-26.py
```

The script will:

1. Simulate coral embryo viability experiments.
2. Estimate statistical power for multiple effect sizes.
3. Determine the minimum number of tubes needed per treatment group.
4. Generate power curve plots.
5. Print summary tables of required replication.

---

## User-Defined Settings

Key parameters can be modified directly in the script:

```python
min_coral_per_tube = 10
max_coral_per_tube = 10

control_p = 0.90

effect_type = "absolute"
effect_sizes = (0.10, 0.20, 0.30)

tube_sd = 0.5

alpha = 0.05
target_power = 0.80
n_sim = 1000
```

### Parameter Descriptions

| Parameter           | Description                                     |
| ------------------- | ----------------------------------------------- |
| control_p           | Baseline viability probability in controls      |
| effect_sizes        | Magnitude of treatment effect to evaluate       |
| effect_type         | "absolute" or "relative" reduction              |
| tube_sd             | Standard deviation of tube-level random effects |
| target_power        | Desired statistical power                       |
| alpha               | Significance level                              |
| n_sim               | Number of Monte Carlo simulations               |
| min_tubes/max_tubes | Range of replication levels evaluated           |

---

## Example Questions Addressed

* How many tubes containing 10 embryos are needed to detect a 10 percentage-point reduction in viability?
* How does tube-to-tube variability affect replication requirements?
* What sample size is needed to achieve 80% statistical power?
* How does detectable effect size influence replication needs?

---

## Repository Structure

```
.
├── PowerAnalysis-CountData-5-27-26.py
├── README.md
└── figures/
```

---

## References

* Quinn, G. P., & Keough, M. J. (2002). Experimental Design and Data Analysis for Biologists.
* Zar, J. H. (2010). Biostatistical Analysis.
* Gelman, A., et al. (2013). Bayesian Data Analysis.
* Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.

---

## Artificial Intelligence Acknowledgment

Artificial intelligence tools, including ChatGPT and Google Gemini, were used to assist in drafting portions of the Python simulation code, reviewing the computational framework, evaluating statistical approaches, and reviewing preliminary versions of the associated Methods text. All AI-generated outputs were critically reviewed, edited, validated, and approved by the author prior to inclusion in this repository.

---

## NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA) or the United States Department of Commerce. All NOAA GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law.

Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply endorsement, recommendation, or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

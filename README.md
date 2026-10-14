# ctDNA--ABM# Agent-Based Model of ctDNA Dynamics

A simple individual-cell agent-based model (ABM), written in R, that simulates how a tumour made of drug-sensitive and drug-resistant cells responds to pulsed treatment, and how this is reflected in circulating tumour DNA (ctDNA) levels.

## What the model does

- Starts with 1,000 tumour cells, 5% of them resistant.
- Each day, every living cell independently and stochastically either undergoes **apoptosis**, **necrosis**, **division**, or nothing.
- **Treatment** is given on days 30, 60 and 90 for 3 days each. During treatment, death probabilities rise sharply for sensitive cells and only slightly for resistant cells.
- **ctDNA** enters the blood from living cells (active release) and from dying cells (apoptosis and, more strongly, necrosis). 80% of the ctDNA pool is cleared each day.
- Runs for 120 days with a fixed random seed (`set.seed(123)`) so results are reproducible.

## Requirements

R (4.1 or newer recommended) and these packages:

```r
install.packages(c("dplyr", "ggplot2", "tidyr", "scales"))
```

## How to run

From the repository folder:

```bash
Rscript ABM_CODE.R
```

or open `ABM_CODE.R` in RStudio, set the working directory to this folder, and click **Source**.

## Outputs

Created inside the folder the script runs in:

| File | Contents |
|---|---|
| `output_data/Simple_ABM_results.csv` | Daily counts of sensitive/resistant cells, deaths, ctDNA level and treatment status |
| `plots/Simple_ABM_total_cells.png` | Total living tumour cells over time |
| `plots/Simple_ABM_sensitive_resistant.png` | Sensitive vs resistant cells |
| `plots/Simple_ABM_ctDNA.png` | Simulated ctDNA dynamics |
| `plots/Simple_ABM_resistant_fraction.png` | Resistant fraction over time |

## Changing the scenario

All parameters sit at the top of the script in clearly labelled sections:

1. **General settings** – simulation length, starting cell count, resistant fraction, treatment days and duration
2. **Biological parameters** – daily division, apoptosis and necrosis probabilities per cell type
3. **Treatment effect** – multipliers applied to death probabilities during treatment
4. **ctDNA release** – release per living or dying cell, and daily clearance fraction

## Notes

Parameter values are illustrative, not fitted to patient data. Simulated ctDNA is in arbitrary units.

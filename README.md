# V2G energy-cost savings

How much can a charging station save on electricity when a driver may cancel all scheduled discharge and keep the promised charging? This repository compares a fixed bidirectional schedule with optimized charging alone under a historical time-of-use tariff. Numerical inputs come from recorded workplace sessions. Battery parameters and bidirectional participation are modeled.

## Run the study

Install the scientific packages with Python 3.13, open [experiments.ipynb](experiments.ipynb) in Jupyter, and select **Run All**. Each run writes a new directory under `runs/` containing the recomputed tables and five figures. All energy-bill experiments use SciPy's HiGHS solvers.

    python -m pip install -r requirements.txt

To reproduce just the figures from the supplied results, run these commands in the repository root. The June 5 figure also solves the two compared schedules.

    python plots/plot_energy_results.py --output-dir figures
    python plots/plot_contract.py --output figures/contract_meter.pdf

The five figures explain the contract, compare bills across 36 fleets, show imports on June 5, show how many withdrawals each fleet can withstand, and report daily savings for July. The [import-limit table](results/common_cap_summary.csv) gives the separate comparison under five shared connection ratings.

## Expected results

| Setting | Charging-only bill | Protected bill | Guaranteed saving |
| --- | ---: | ---: | ---: |
| 36 fleets, one withdrawal per fleet-day, unrestricted import | $1,436.96 | $1,392.90 | $44.06 (3.07%) |
| Same fleets, minimum feasible charging-only import limit | $1,570.66 | $1,568.96 | $1.69 (0.11%) |
| All 31 days of July, one withdrawal per day, unrestricted import | $1,416.71 | $1,328.79 | $87.92 (6.21%) |

Values are rounded independently. Each protected bill is the highest bill among admissible withdrawals at the optimized timetable. The July cohort overlaps with the main study. The tariff excludes demand charges and gives zero credit for exports.

## Command-line replication

From the repository root, the following commands recompute the numerical results in the supplied output folders. The notebook runs the same experiments in a separate folder. The tested environment was Python 3.13.2, NumPy 2.2.6, SciPy 1.17.1 and Matplotlib 3.10.0.

    python run_study.py all
    python src/validate_threshold.py --fleets inputs/study --out results
    python run_cap_sensitivity.py
    python mode_replication/mode_check.py --fleets inputs/study --published results/q_sweep.csv --out mode_replication/mode_checks.csv --limit 20
    python mode_replication/surplus_check.py --out mode_replication/surplus_mode_check.csv
    python payment_replication/discharge_allowance.py results --out payment_replication/discharge_allowance.csv

The threshold sweep covers 540 fleet/budget pairs. The [mode checks](mode_replication/README.md) compare continuous and binary-mode schedules. The [common-limit experiment](README_caps.md) solves both policies under each of 180 matched fleet/limit combinations. To independently rebuild the 67 numerical inputs from the pinned processed Caltech export, see [reconstruction](reconstruction/README.md).

## Inputs

The numerical inputs represent fleets constructed from a processed export of Caltech workplace sessions in ACN-Data. Observed fields are connection time, disconnection time and delivered grid-side energy. A session is eligible if it is distinct, is entirely within the local calendar day, delivered at least 2 kWh and can deliver the recorded energy in its available window at the assigned charging power. Local times are interpreted in America/Los_Angeles. The time grid splits at hour boundaries, connections and disconnections. Session and registry identifiers are omitted. An output vehicle index is local to its fleet and begins at one.

**inputs/study** contains 36 fleet-days (504 sessions) from April through September 2019. Each month contributes six distinct days with 6, 12, 24, 6, 12 and 24 vehicles. The construction chose eligible dates in order and selected sessions by a fixed hash (selection seed 20261114). **inputs/july** contains all eligible sessions on each of the 31 days of July 2019. Of 811 arrivals in the processed export, 158 failed the energy screen, 39 extended outside the day and 36 could not receive the recorded energy at the assigned charging power. The 578 retained sessions overlap with the July subsample in **inputs/study**. **inputs_sha256.csv** records a SHA-256 digest for every numerical input.

In each JSON file, availability is a vehicle-by-interval array, duration is in hours, and boundaries_seconds gives seconds since local midnight. Pc and Pd are grid-side charging and discharge limits in kW. E0, Emin, Emax and Etarget are battery energies in kWh. eta_c and eta_d are dimensionless efficiencies. L is a constant 10 kW background load. Both scaling arrays, Ulo and Uhi, equal one.

Battery parameters and permissions are modeled. Initial energies are 15, 25 or 35 kWh and lower bounds are zero. The target equals initial energy plus charging efficiency times recorded grid energy. Upper battery limits add 0, 0.5, 5 or 15 kWh above the target. Charging ratings are 3.6, 7.2 or 11 kW, discharge ratings are 2, 4 or 7.2 kW, and charging and discharging efficiencies are drawn from (0.90, 0.93, 0.96) and (0.88, 0.92, 0.96). The numerical arrays preserve the fixed assignment to each session. The observed data contain no measurements of V2G decisions or withdrawals in this study.

The original session data are available from [Caltech ACN-Data](https://ev.caltech.edu/dataset). Please cite Z. J. Lee, T. Li and S. H. Low, [ACN-Data: Analysis and Applications of an Open EV Charging Dataset](https://doi.org/10.1145/3307772.3328313), ACM e-Energy 2019, when using these data. The 36 sampled fleet inputs follow the temporal cohort in the [capacity study repository](https://github.com/dhducpro/v2g-capacity).

## Prices

The energy prices follow the historical Southern California Edison TOU-EV-4 benchmark reported by Lee, Pang and Low in their Table I and simulation description. From June through September the off-, mid- and on-peak rates are 0.0562, 0.0925 and 0.267 USD/kWh. The other months use 0.061, 0.075 and 0.087 USD/kWh. Weekday on-peak hours are 12:00–18:00 local time and mid-peak hours are 08:00–12:00 and 18:00–23:00. All other hours and entire weekends are off-peak. Weekday holidays follow the weekday rule. The prices provide a benchmark for purchased energy. The campus's actual retail contract is unknown. [Pricing EV Charging Service with Demand Charge](https://pscc-central.epfl.ch/repo/papers/2020/400.pdf).

## Interpreting the results

| Output in results/ | Contents |
| --- | --- |
| daily_monetary_by_fleet.csv | One row per study fleet. The prefixes charge_only, nominal_first and robust_first identify the independently optimized policies. Each includes nominal and worst energy cost and physical audit residuals. Costs include background demand. |
| same_fleets_by_season.csv | Each of the 36 fixed fleets evaluated under both seasonal tariffs. Its actual weekday/weekend status is preserved. |
| july_by_day.csv, july_totals.csv | Daily optimized charging-only and robust costs, and monthly sums. One withdrawal is allowed on each date. The energy objective and constraints separate across dates. |
| common_cap_summary.csv, common_cap_by_fleet.csv | Pooled and fleet-level comparisons under the same import limit for both policies. |
| driver_exposures.csv | Dimensionless individual-vehicle scores from the criterion's linear program. |
| fleet_thresholds.csv | Predicted and observed first withdrawal budgets at which guaranteed savings vanish. |
| q_sweep.csv | Clipped-meter cost bounds, realized savings, feasibility residuals and independent checks for all 540 fleet/budget pairs. A blank explicit-check cell means this point was solved by separation without the additional fully enumerated comparison. |

For the one-withdrawal dispatch runs, ties in the primary energy objective are resolved by minimizing the other energy objective followed by the protected peak. A tolerance of 1e-8 USD accounts for rounding between optimization stages. **src/energy.py** checks the recovered physical schedules and recomputes the bills for each admissible withdrawal. **src/validate_threshold.py** compares lower bounds from a linear program against the highest attainable clipped bills and adds violated constraints as needed. It also solves fully enumerated scenario models at 132 comparison points. Savings above 1e-6 USD are counted as positive.

Three checks on the supplied outputs are a guaranteed saving of 44.06 USD across the 36 study fleets at one withdrawal, 87.92 USD over the 31 days of July, and agreement on the sign of savings at all 540 fleet/budget pairs. Floating-point values and schedules among equivalent optima may vary slightly between solver versions.

## Original model controls

The earlier model combined energy trading with throughput and quadratic smoothing. Its [original-model control archive](https://github.com/dhducpro/v2g-capacity/blob/2f7c7586ac9f7852be1657336f0ab4339b1add0a/reviewer-control-records.zip) contains the code and run-level records for that separate objective. The energy-procurement study in this repository measures costs in US dollars.

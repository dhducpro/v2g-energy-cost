# Independent explicit-mode verification of the revised energy model

These files verify the uncapped tariff-priced energy objective with one protected whole-session withdrawal. The inputs and q_sweep.csv belong to the accompanying main replication directory. Each case is solved twice with an explicit scenario import meter, once as a linear program and once with one binary charging/discharging variable per connected vehicle and interval. Solver-independent audits recompute purchased energy from the sign of realized net demand.

From this directory, pass the path to the extracted replication inputs and reference CSV. For example, if this directory is inside the main replication directory:

    python mode_check.py --fleets ../inputs/study --published ../results/q_sweep.csv --out mode_checks.csv --limit 20
    python surplus_check.py --out surplus_mode_check.csv

The first command overwrites mode_checks.csv with results and measured wall time on the executing machine. The independent verifier uses its own sparse matrix, interval-end energy state variables, a common signed station meter in each interval, an import epigraph for each complete withdrawal scenario, and a tariff constructed directly from the study day. No capacity or nominal-bill tie-breaking is used. Every solver was run with one HiGHS thread. The requested MILP relative gap was zero, with mixed-integer feasibility tolerance 1e-8 and LP primal feasibility tolerance 1e-9. The audit tolerance on the dollar objective and the binary schedule's battery and overlap residuals was 1e-6. Measurements exclude matrix building and auditing. Python 3.13.2, NumPy 2.2.6, SciPy 1.17.1 and Windows 11 were used. The 20-s time limit applies to each solve separately. SciPy status 0 denotes optimal and status 1 denotes a limit.

The surplus control has two 1-hour intervals, net background demand -4 kW then 3 kW, positive purchase prices 0.10 then 0.30 USD/kWh, and zero export credit. The two EVs both start and must depart at 2 kWh, with bounds [0, 2.9] kWh, powers up to 1 kW and charging/discharging efficiencies 0.9/0.9. Every vehicle is connected in both intervals. The surplus case deliberately lies outside the positive-background assumption of the score theorem. Its LP and explicit-mode models should agree on 0.657 USD under one withdrawal, against 0.900 USD from optimized charging only.

Reference outputs are included for both commands. CPU times vary by machine. MILP full-delivery PV export and minimum signed net demand describe one optimal returned schedule and need not characterize all optimal timetables.

# Matched site-capacity sensitivity

From the extracted replication root, run

    python run_cap_sensitivity.py

The script uses only `inputs/study` and `src/energy.py`. It first solves an independent linear program for each day's minimum charging-only peak station import. It then solves a new charging-only schedule and a new protected bidirectional schedule at each of five common physical import limits, with one whole-session withdrawal per day. The minimum includes the background load. The four finite limits equal 1, 1.10, 1.25 and 1.50 times that day's minimum charging-only peak. Unit net power coefficients apply to all 36 study inputs.

Outputs are `results/common_cap_by_fleet.csv` (180 paired observations) and `results/common_cap_summary.csv` (five rows). Percentage savings divide aggregate savings by the charging-only bill under the same limit. A monetary saving counts as positive above USD 1e-6. The schedule solver audits battery feasibility, achieved peak import and each scenario's bill. When the limit is relaxed, both optimum bills must weakly decrease, but their difference need not increase for each fleet.

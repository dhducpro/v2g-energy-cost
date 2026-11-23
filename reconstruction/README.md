# Rebuild the numerical inputs from the processed session records

The archive contains a deterministic reconstruction from the processed Caltech session export to the vehicle inputs used in the 36-fleet study and the full July energy-cost experiment. The source export is the historical nine-column, 27,558-record snapshot. Only Start, End and Energy are used here. The original ACN session database is described by [Lee, Li and Low (2019)](https://doi.org/10.1145/3307772.3328313) and available from [Caltech](https://ev.caltech.edu/dataset). This utility starts from the pinned processed export. It does not claim to recreate that nine-column snapshot from a fresh API query.

Download the original [study-data.zip](https://github.com/dhducpro/v2g-capacity/raw/aebe6f8a49243be42880fb0471a46627807c7a31/study-data.zip) from the authors' [capacity-study repository](https://github.com/dhducpro/v2g-capacity/tree/aebe6f8a49243be42880fb0471a46627807c7a31). Its SHA-256 digest is `484b607820ca818757eef6b9ecc7c09b8d68149a0cb4123e66100ae23f7c0a3c`. The script checks both the archive and each input member before use. It reads the three files below in memory.

- `replication/cohort_reconstruction/inputs/caltech.csv` (processed sessions)
- `replication/cohort_reconstruction/inputs/frozen_selected_sessions.csv` (144 earlier selections)
- `replication/cohort_reconstruction/inputs/original_excluded_sessions.csv` (18 earlier selections)

The exclusion sets preserve the original study protocol. All 162 dates of connection for these sessions fall in January--March; they add no July exclusions. Session keys built from UTC connection and departure timestamps and reported grid energy are used internally to identify physical duplicates and apply deterministic sampling and parameter assignments. Emitted JSONs carry numerical input arrays and the horizon date only. No upstream session identifier or individual registry is emitted. The processed export and earlier manifests are not redistributed with this energy-cost archive.

From the replication directory, run

    python reconstruction/rebuild_inputs.py --archive /path/to/study-data.zip

An optional `--out rebuilt` writes a second copy under `rebuilt/inputs/study` and `rebuilt/inputs/july` after all checks pass. Existing files are preserved. The code uses the Python standard library and an IANA time-zone database for America/Los_Angeles (install `tzdata` if the system lacks it).

The reconstruction checks all 67 SHA-256 digests against `inputs_sha256.csv`. They cover the complete numerical files and their sanitized metadata. The 36 sampled fleet-days contain 504 selected sessions. All 31 July fleet-days contain 578 eligible sessions out of 811 recorded arrivals. The July eligibility checks remove 158 sessions at the energy screen, 39 that cross midnight and 36 that cannot finish at their assigned charging rate. The complete July set overlaps the July portion of the 36-fleet sample.

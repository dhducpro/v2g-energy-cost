# V2G energy-cost savings

This study asks whether a charging station can still save on purchased electricity when drivers cancel discharge but keep every scheduled charge. It uses observed Caltech charging sessions and a historical time-of-use tariff.

## Reproduce the results

From the repository root, run these commands with Python 3.13.

```sh
python -m pip install -r requirements.txt jupyterlab
python -m jupyterlab experiments.ipynb
```

Select **Run All**. The notebook checks the 67 input files and writes fresh tables and five figures under `runs/`. Across 36 fleets, allowing one cancellation per day gives a guaranteed energy-cost saving of $44.06 (3.07%) with unrestricted imports.

[Reference results](results/) · [Data source and reconstruction](reconstruction/README.md) · [Import limits](README_caps.md) · [Battery-mode checks](mode_replication/README.md) · [Original-model controls](https://github.com/dhducpro/v2g-capacity/blob/2f7c7586ac9f7852be1657336f0ab4339b1add0a/reviewer-control-records.zip)

#!/usr/bin/env python3

import numpy as np

data = np.loadtxt("checkReconstructedTmunu.dat", ndmin=2)
if data.size == 0 or data.shape[1] != 3 or not np.isfinite(data).all():
    raise SystemExit("Invalid or nonfinite reconstruction data")

mean_residual = np.mean(np.abs(data[:, 2]))
if not np.isfinite(mean_residual) or mean_residual > 0.001:
    raise SystemExit(f"Diff: {mean_residual}")

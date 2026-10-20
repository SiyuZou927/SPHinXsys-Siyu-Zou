"""Run with ParaView pvpython; verify the replica before integrating fluid motion."""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import vtk_to_numpy

new_csv, reference, result = map(Path, sys.argv[1:4])
geometry = pd.read_csv(new_csv)
new = geometry[["initial_x", "initial_y", "initial_z"]].to_numpy()
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(str(reference))
reader.Update()
old = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
report = {"new_count": len(new), "reference_count": len(old),
          "reference": str(reference), "new_geometry": str(new_csv),
          "shoulder_mask_count": int(geometry.annular_shoulder_mask.sum()),
          "station_range_m": [float(geometry.station_m.min()), float(geometry.station_m.max())],
          "axial_bin_counts": geometry.bin.value_counts().sort_index().to_dict()}
report["constrained_count"] = int(geometry.simbody_constrained.sum())
report["volume_range_m3"] = [float(geometry.volume_m3.min()), float(geometry.volume_m3.max())]
report["constrained_mass_assuming_rho_2_063"] = float(
    geometry.loc[geometry.simbody_constrained == 1, "volume_m3"].sum()*2.063)
report["volume_verified_against_old_vtp"] = False
if new.shape == old.shape:
    # Reload and VTP retain material indices; checking indexwise is stronger than
    # a nearest-neighbor fit and catches a changed particle realization.
    report["max_indexwise_distance_m"] = float(np.linalg.norm(new-old, axis=1).max())
    report["mean_indexwise_distance_m"] = float(np.linalg.norm(new-old, axis=1).mean())
    report["passed"] = report["max_indexwise_distance_m"] < 2e-6
else:
    report["passed"] = False
result.write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
raise SystemExit(0 if report["passed"] else 2)

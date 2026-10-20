"""Read-only cross-check of archived VTP poses and Summary body states."""
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from run_diagnostic import ROOT, REFERENCE


def positions(path):
    root = ET.parse(path).getroot()
    return np.fromstring(root.find("./PolyData/Piece/Points/DataArray").text, sep=" ").reshape(-1, 3)


p0 = positions(REFERENCE / "FloatingObject_0000000000.vtp")
c0 = p0.mean(axis=0)
x = p0-c0
b = pd.read_csv(REFERENCE / "SummaryOutput3D.dat", sep=r"\s+")
rows = []
for path in sorted(REFERENCE.glob("FloatingObject_*.vtp")):
    label = int(path.stem.rsplit("_", 1)[1])
    if label == 0: continue
    time = round(label/1000)*.001
    if time > .008: continue
    points = positions(path)
    c = points.mean(axis=0)
    y = points-c
    u, _, vt = np.linalg.svd(x.T@y)
    rotation = vt.T@u.T
    translation = c-rotation@c0
    fit = p0@rotation.T+translation
    center = np.array([np.interp(time, b.iloc[:, 0], b[col]) for col in
                       ["MassCenter_X", "MassCenter_Y", "MassCenter_Z"]])
    angle = np.interp(time, b.iloc[:, 0], b.RotationAngle_Z)
    inferred = rotation.T@(center-translation)
    rows.append(dict(time_s=time, vtp_label=label,
                     rigid_fit_max_m=float(np.linalg.norm(fit-points, axis=1).max()),
                     inferred_reference_center_x=float(inferred[0]),
                     inferred_reference_center_y=float(inferred[1]),
                     inferred_reference_center_z=float(inferred[2]),
                     inferred_vs_initial_cloud_center_m=float(np.linalg.norm(inferred-c0)),
                     vtp_cloud_vs_summary_center_m=float(np.linalg.norm(c-center)),
                     rigid_rotation_rad=float(np.arctan2(rotation[1, 0], rotation[0, 0])),
                     summary_rotation_change_rad=float(angle-np.deg2rad(-20))))
out = ROOT / "legacy_pose_crosscheck.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows).to_string(index=False))

force = b.PressureForceGlobal_X+b.ViscousForceGlobal_X
unloaded = b[force.eq(0) & b.HydrodynamicMoment_Z.eq(0)].copy()
unloaded = unloaded[unloaded.iloc[:, 0] < .001]
vx0 = 70*np.cos(np.deg2rad(20))
offset = (unloaded.MassCenter_X-c0[0])/vx0-unloaded.iloc[:, 0]
report = dict(zero_force_early_records=len(unloaded),
              inferred_clock_offset_s=offset.tolist(),
              median_clock_offset_s=float(offset.median()),
              one_initial_acoustic_dt_estimate_s=.6*1.3*.005/700,
              explanation="A clock offset is a diagnostic possibility; this does not prove its source.")
(ROOT / "legacy_clock_crosscheck.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))

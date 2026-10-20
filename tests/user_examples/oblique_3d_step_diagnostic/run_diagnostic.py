"""One isolated baseline run; production inputs and manuscript remain untouched."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

import numpy as np
import pandas as pd

SOURCE = Path(__file__).resolve().parent
WORKSPACE = SOURCE.parents[3]
BUILD = WORKSPACE / "build/tests/user_examples/oblique_3d_step_diagnostic/RelWithDebInfo"
ROOT = Path(r"G:\SPHdata\oblique_water_entry\diagnostic_3d_step_20261018")
REFERENCE = Path(r"G:\SPHdata\oblique_water_entry\3D C=700\gama25\5MM\output")
PVPYTHON = Path(r"D:\Program Files\ParaView 6.0.1\bin\pvpython.exe")
RUN = ROOT / "baseline_v70_c700_gamma25_5mm"
PREFLIGHT = ROOT / "initialization_preflight"


def sha(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, obj):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    temp.replace(path)


def state(stage, **kwargs):
    save(ROOT / "run_status.json", dict(stage=stage, updated_utc=datetime.now(timezone.utc).isoformat(),
                                       runner_pid=os.getpid(), **kwargs))


def observable_onset(t, signed_force):
    compression = np.maximum(-signed_force, 0) / (.5*70.0**2*np.pi*(.0211/2)**2)
    crossing = np.flatnonzero((compression[:-1] >= .001) & (compression[1:] >= .001))
    if not len(crossing):
        return {"time_s": None, "left_censored": False, "detected": False}
    i = int(crossing[0])
    if i == 0:
        onset = float(t[0])
    else:
        onset = float(t[i-1]+(.001-compression[i-1]) /
                      (compression[i]-compression[i-1])*(t[i]-t[i-1]))
    return {"time_s": onset, "left_censored": i == 0, "detected": True}


def prepare(reload_file):
    exe = BUILD / "oblique_3d_step_diagnostic.exe"
    if not exe.exists():
        raise FileNotFoundError(exe)
    if ROOT.exists() and set(p.name for p in ROOT.iterdir()) != {"reconstructed_input"}:
        raise FileExistsError(f"Do not overwrite an existing diagnostic package: {ROOT}")
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT / "source_snapshot").mkdir()
    for path in SOURCE.iterdir():
        if path.is_file():
            shutil.copy2(path, ROOT / "source_snapshot" / path.name)
    provenance = []
    for directory in [RUN, PREFLIGHT]:
        for name in ["input", "reload", "bin", "logs"]:
            (directory / name).mkdir(parents=True)
        for path in BUILD.iterdir():
            if path.suffix.lower() in [".exe", ".dll"]:
                shutil.copy2(path, directory / "bin" / path.name)
        # Dependencies already supplied by the original successfully built case.
        dependencies = WORKSPACE / "build/tests/3d_examples/test_3d_zsy_oblique_water_entry_exit/bin/RelWithDebInfo"
        for path in dependencies.glob("*.dll"):
            dest = directory / "bin" / path.name
            if not dest.exists():
                shutil.copy2(path, dest)
        stl = Path(r"E:\zsy study\oblique_water_entry\3Dmodel\M1TR4Z02H05.stl")
        shutil.copy2(stl, directory / "input" / stl.name)
        shutil.copy2(reload_file, directory / "reload/FloatingObject_rld.xml")
        for path in directory.rglob("*"):
            if path.is_file():
                provenance.append(dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size))
    config = dict(V0_mps=70.0, c0_mps=700.0, gamma=25.0, dp_fluid_m=0.005,
                  solid_adaptation_ratio=2.0, trajectory_angle_deg=-20.0,
                  pitch_angle_deg=-20.0, end_time_s=0.008,
                  force_interval_s=0.00004, vtp_interval_s=0.001,
                  reload_source=str(reload_file), reference=str(REFERENCE),
                  production_source_edited=False, manuscript_edited=False,
                  diagnostic_only=True, initial_state_reconstructed_from_vtp=True,
                  uniform_solid_volume_assumed_m3=.0025**3,
                  author_reconstruction_approved_date="2026-10-18")
    save(ROOT / "configuration.json", config)
    save(ROOT / "input_sha256.json", provenance)
    state("prepared")


def preflight():
    env = os.environ.copy()
    env["SPH_DIAGNOSTIC_INITIALIZE_ONLY"] = "1"
    state("initializing_preflight")
    with (PREFLIGHT / "logs/solver.log").open("w", encoding="utf-8") as log:
        code = subprocess.run([str(PREFLIGHT / "bin/oblique_3d_step_diagnostic.exe")],
                              cwd=PREFLIGHT, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    if code:
        raise RuntimeError(f"Preflight exit code {code}")
    validate_preflight_data()


def validate_preflight_data():
    subprocess.run([str(PVPYTHON), str(SOURCE / "check_initial_geometry.py"),
                    str(PREFLIGHT / "output/DiagnosticGeometry.csv"),
                    str(REFERENCE / "FloatingObject_0000000000.vtp"),
                    str(ROOT / "initial_geometry_validation.json")], check=True)
    subprocess.run([sys.executable, str(SOURCE / "check_legacy_kinematics.py")], check=True)
    subprocess.run([sys.executable, str(SOURCE / "validate_mass_motion.py"), str(ROOT)], check=True)
    state("preflight_passed")


def analyze():
    d = pd.read_csv(RUN / "output/PressureDecomposition.csv")
    s = pd.read_csv(RUN / "output/SummaryOutput3D.dat", sep=r"\s+")
    w = pd.read_csv(RUN / "output/AppliedWrench.csv")
    if not len(s) or abs(float(s.iloc[-1, 0]) - 0.008) > 1e-9:
        raise ValueError("Simulation has not reached the requested 8 ms endpoint")
    expected = 200
    if len(s) != expected or len(w) != expected:
        raise ValueError(f"Expected {expected} output events, got {len(s)}, {len(w)}")
    for name, table in [("summary", s), ("applied_wrench", w), ("decomposition", d)]:
        if not np.isfinite(table.select_dtypes(include=[np.number]).to_numpy()).all():
            raise ValueError(f"Non-finite diagnostic values in {name}")
    if not np.all(np.diff(s.iloc[:, 0]) > 0) or not np.all(np.diff(w.time_s) > 0):
        raise ValueError("Output time must be strictly increasing")
    counts = d.groupby(["stage", "region"]).size()
    if len(counts) != 28 or not (counts == expected).all():
        raise ValueError("Incomplete region/stage diagnostic history")
    all_force = d[d.region == "all"]
    scale = max(1.0, float(all_force.pressure_axial.abs().max()))
    max_closure = float(d.closure_norm.max())
    partition_errors = []
    for (stage_name, t), g in d.groupby(["stage", "time_s"]):
        total = g[g.region == "all"].iloc[0]
        bins = g[g.region.str.startswith("axial_bin_")]
        for col in ["pressure_axial", "central_axial", "acceleration_axial", "riemann_axial", "viscous_axial"]:
            partition_errors.append(abs(float(bins[col].sum() - total[col])))
    b = pd.read_csv(REFERENCE / "SummaryOutput3D.dat", sep=r"\s+")
    t = s.iloc[:, 0].to_numpy()
    baseline = np.interp(t, b.iloc[:, 0], b.TotalForceLocal_X)
    replica = s.TotalForceLocal_X.to_numpy()
    comparison = pd.DataFrame(dict(time_s=t, baseline_axial=baseline, diagnostic_axial=replica,
                                   residual_axial=replica-baseline))
    comparison.to_csv(ROOT / "baseline_comparison.csv", index=False)
    fresh = d[d.stage == "fresh_output"]
    rows = []
    columns = ["pressure_axial", "central_axial", "acceleration_axial", "riemann_axial", "viscous_axial", "total_axial"]
    for region, g in fresh.groupby("region"):
        before = g[g.time_s.between(.0050, .0056)]
        after = g[g.time_s.between(.0064, .0070)]
        for col in columns:
            rows.append(dict(region=region, term=col, before_mean=float(before[col].mean()),
                             after_mean=float(after[col].mean()),
                             change=float(after[col].mean()-before[col].mean())))
    pd.DataFrame(rows).to_csv(ROOT / "step_window_decomposition.csv", index=False)
    report = dict(summary_records=len(s), applied_wrench_records=len(w),
                  decomposition_records=len(d), end_time_s=float(t[-1]),
                  max_decomposition_closure_norm=max_closure,
                  max_partition_sum_error=float(max(partition_errors)),
                  baseline_axial_correlation=float(np.corrcoef(replica, baseline)[0, 1]),
                  baseline_axial_rmse=float(np.sqrt(np.mean((replica-baseline)**2))),
                  baseline_peak_relative_rmse=float(np.sqrt(np.mean((replica-baseline)**2))/max(abs(baseline))),
                  diagnostic_load_onset=observable_onset(t, replica),
                  baseline_load_onset=observable_onset(b.iloc[:, 0].to_numpy(), b.TotalForceLocal_X.to_numpy()),
                  numerical_closure_passed=max_closure < 1e-10*scale and max(partition_errors) < 1e-10*scale,
                  interpretation="Term and region differences are diagnostics, not a unique causal identification.")
    save(ROOT / "completion_checks.json", report)
    if not report["numerical_closure_passed"]:
        raise ValueError("Diagnostic decomposition failed numerical closure")
    return report


def run():
    validation = json.loads((ROOT / "initial_geometry_validation.json").read_text())
    if not validation["passed"]:
        raise ValueError("Initial geometry verification is required")
    mass_motion = json.loads((ROOT / "initial_mass_motion_validation.json").read_text())
    if not mass_motion["passed"]:
        raise ValueError("Mass and motion consistency verification is required")
    if (RUN / "output/SummaryOutput3D.dat").exists():
        raise FileExistsError("Never overwrite an already-started run")
    state("running", end_time_s=.008)
    env = os.environ.copy()
    env.pop("SPH_DIAGNOSTIC_INITIALIZE_ONLY", None)
    with (RUN / "logs/solver.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen([str(RUN / "bin/oblique_3d_step_diagnostic.exe")],
                                   cwd=RUN, env=env, stdout=log, stderr=subprocess.STDOUT)
        state("running", solver_pid=process.pid, end_time_s=.008)
        code = process.wait()
    if code:
        raise RuntimeError(f"Solver exit code {code}")
    state("analyzing")
    report = analyze()
    state("complete", checks=report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["prepare", "preflight", "validate", "run", "analyze"])
    parser.add_argument("--reload", type=Path)
    args = parser.parse_args()
    if args.action == "prepare": prepare(args.reload)
    elif args.action == "preflight": preflight()
    elif args.action == "validate": validate_preflight_data()
    elif args.action == "run": run()
    else: print(json.dumps(analyze(), indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        if ROOT.exists(): state("failed", traceback=traceback.format_exc())
        raise

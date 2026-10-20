"""Synthetic checks for postflight acceptance, not simulation evidence."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
import run_diagnostic as run


class PostflightTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        run.ROOT = self.root
        run.RUN = self.root / "run"
        run.REFERENCE = self.root / "reference"
        (run.RUN / "output").mkdir(parents=True)
        run.REFERENCE.mkdir()
        t = np.arange(1, 201)*.00004
        summary = pd.DataFrame({"Time[s]": t, "TotalForceLocal_X": -t*120})
        summary.to_csv(run.RUN / "output/SummaryOutput3D.dat", sep=" ", index=False)
        summary.to_csv(run.REFERENCE / "SummaryOutput3D.dat", sep=" ", index=False)
        pd.DataFrame({"time_s": t, "dt_s": np.full(200, 1e-6)}).to_csv(
            run.RUN / "output/AppliedWrench.csv", index=False)
        rows = []
        for stage in ["cached_before_density", "fresh_output"]:
            for time in t:
                for region in [f"axial_bin_{i}" for i in range(12)]+["all", "annular_shoulder_mask"]:
                    force = -time*120 if region in ["all", "axial_bin_0"] else 0.0
                    rows.append(dict(stage=stage, time_s=time, region=region, pressure_axial=force,
                                     central_axial=force*.7, acceleration_axial=force*.1,
                                     riemann_axial=force*.2, viscous_axial=0.0,
                                     total_axial=force, closure_norm=0.0))
        self.path = run.RUN / "output/PressureDecomposition.csv"
        pd.DataFrame(rows).to_csv(self.path, index=False)

    def tearDown(self):
        self.temp.cleanup()

    def test_complete(self):
        report = run.analyze()
        self.assertTrue(report["numerical_closure_passed"])
        self.assertAlmostEqual(report["baseline_axial_correlation"], 1.0)

    def test_bad_closure_rejected(self):
        d = pd.read_csv(self.path)
        d.loc[0, "closure_norm"] = 1.0
        d.to_csv(self.path, index=False)
        with self.assertRaises(ValueError): run.analyze()

    def test_missing_record_rejected(self):
        d = pd.read_csv(self.path).iloc[:-1]
        d.to_csv(self.path, index=False)
        with self.assertRaises(ValueError): run.analyze()


if __name__ == "__main__":
    unittest.main()

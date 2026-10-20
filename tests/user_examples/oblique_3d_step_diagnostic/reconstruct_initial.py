"""Rebuild a separate reload from archived initial VTP; never relax or overwrite."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np

from run_diagnostic import ROOT, REFERENCE

source = REFERENCE / "FloatingObject_0000000000.vtp"
destination = ROOT / "reconstructed_input"
if destination.exists():
    raise FileExistsError(destination)
document = ET.parse(source)
piece = document.find("./PolyData/Piece")
points_array = piece.find("./Points/DataArray")
assert points_array.attrib.get("Format", "ascii").lower() == "ascii"
points = np.fromstring(points_array.text, sep=" ").reshape(-1, 3)
ids_array = piece.find("./PointData/DataArray[@Name='OriginalID']")
ids = np.fromstring(ids_array.text, sep=" ", dtype=np.int64)
assert len(points) == len(ids) == int(piece.attrib["NumberOfPoints"]) == 57399
assert np.isfinite(points).all() and len(np.unique(ids)) == len(ids)
volume = .0025**3
root = ET.Element("particles")
for original_id, point in zip(ids, points):
    ET.SubElement(root, "particle", OriginalID=str(int(original_id)),
                  VolumetricMeasure=format(volume, ".17g"),
                  Position=", ".join(format(float(v), ".17g") for v in point))
destination.mkdir(parents=True)
output = destination / "FloatingObject_rld.xml"
ET.indent(root)
ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
back = ET.parse(output).getroot()
reloaded = np.array([np.fromstring(p.attrib["Position"], sep=",") for p in back])
assert np.array_equal(reloaded, points)
mass = len(points)*volume*2.063
center = points.mean(axis=0)
d = points-center
inertia = mass*(np.eye(3)*np.mean(np.sum(d*d, axis=1))-d.T@d/len(d))
report = dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              output=str(output), output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
              count=len(points), coordinates="Original ASCII VTP decimals retained in order, no relaxation",
              volume_assumption_m3=volume, volume_not_present_in_source_vtp=True,
              density_model_units=2.063, full_cloud_mass_model_units=mass,
              full_cloud_center_m=center.tolist(), full_cloud_inertia_model_units_m2=inertia.tolist(),
              xml_roundtrip_exact=True, author_approved_date="2026-10-18")
(destination / "reconstruction.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))

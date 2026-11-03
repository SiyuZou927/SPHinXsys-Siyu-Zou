# STL geometry for the 3D oblique water-entry example

## Included geometry

`regression_test_tool/M1TR4Z02H05.stl` is the unmodified binary STL referenced by
`test_3d_water_entry_exit.cpp` as `./input/M1TR4Z02H05.stl`.

The existing case CMake configuration copies files from `regression_test_tool/`
to the executable's `bin/input/` directory when CMake configures the case.
Reconfigure CMake after obtaining this addition. Alternatively, copy this exact
file manually to `input/M1TR4Z02H05.stl` relative to the run working directory.
Do not overwrite unrelated input files.

The source uses zero translation and scale factor1.0 for this STL. Coordinates
are already in the simulation's metre-scale coordinate system and include the
initial position and orientation. Do not convert millimetres to metres or apply
another20-degree rotation on loading. The principal-axis angle measured from the
mesh is approximately -19.98 degrees in the x-y plane; this is a mesh diagnostic,
not a replacement for the nominal -20-degree setup.

| Property | Value |
| --- | --- |
| File size | 22884 bytes |
| Triangles | 456 |
| Minimum x,y,z / m | 0.006705347,0.528069794,0.170767814 |
| Maximum x,y,z / m | 0.467146605,0.718930900,0.229232177 |
| SHA-256 | `75777eaa50ce1ae354ddcfc830dac9d47edd7617fae5a7f213a37149a3e35738` |

## Generate the solid reload before a fresh run

The example defaults to loading a relaxed solid-particle file. The STL alone is
not that particle reload. After building the3D target, run from its `bin/`
working directory so that `./input/` is resolved correctly. With no existing
reload, first generate one using the example's original relaxation path:

```text
test_3d_zsy_oblique_water_entry_exit --relax true --reload false --restart_step 0
```

Then perform a fresh simulation with the generated reload:

```text
test_3d_zsy_oblique_water_entry_exit --relax false --reload true --restart_step 0
```

Use `.exe` on Windows if needed. Do not reuse an unrelated solid reload or a
different geometry's particle file. A fresh randomized relaxation need not
reproduce a historical particle realization bit for bit. No archived reload,
restart, simulation results or regression-output XML is added by this commit.
These instructions document the source workflow; this asset-publication step
did not rebuild or rerun the simulation.

## Important: source defaults differ from the blog case

The blog and [SPHinXsys Discussion1131](https://github.com/Xiangyu-Hu/SPHinXsys/discussions/1131)
describe a historical V70/c700 run with wetting parameter gamma=25. The linked
organized source at commit`b657f4870c1617c3e1e4cffa56b307812ccd9e36`
has different default parameters:

| Setting | Organized-source default | Blog case |
| --- | --- | --- |
| Initial speed | V50:50m/s | V70:70m/s |
| Sound speed | 10*initial_speed=500m/s | 700m/s |
| Wetting diffusion multiplier | 0*dx^2 | 25*dx^2 |
| Fluid particle spacing | 0.005m | 0.005m |
| Nominal solid spacing | 0.0025m from adaptation ratio2 | 0.0025m |
| Nominal entry angle | -20degrees | -20degrees |

Adding the STL closes the missing-geometry dependency; it does not make the
unmodified default source an exact reproduction of the blog's load history.
For a comparison with that historical case, use a separate configuration with
the documented V70/c700/gamma25 values and record all other source and particle
initialization differences. This commit deliberately changes no solver code,
material setting, force definition or time-integration algorithm. It does not
certify that the historical numerical force step is physical.

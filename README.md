# Automatic base-parameter regrouping for robot dynamic equations

MATLAB implementation of the symbolic–numerical procedure that derives **explicit, certified base
(minimal inertial) parameters** directly from the symbolic joint-torque equations of a serial
manipulator.

> M. A. Yusuf, N. M. Alyazidi, Md Shafiullah, A. Nasir,
> *Automatic Symbolic–Numerical Derivation of Explicit Base Parameters from Robot Dynamic Equations*
> (submitted to Mechanism and Machine Theory; citation to be added after publication).
>
> Repository: https://github.com/Muhammad-Adel-Yusuf/base-parameter-regrouping

## What it does

Given torque expressions that are linear in parameter monomials, the procedure

1. extracts the parameter monomials and builds a common regressor,
2. finds the structural rank by SVD at random states,
3. regroups the dependent monomials by column-pivoted QR into explicit base parameters,
4. verifies the result numerically and **certifies it exactly** by polynomial reduction
   modulo `sin² + cos² − 1`.

Link dimensions can stay **symbolic** (the result then holds for any dimension values) or be
substituted with known values. Rotor inertia and friction (`Ia`, `Fv`, `Fc`) are handled as
additional symbols without any change. The method accepts torque equations from any derivation
(Euler–Lagrange or Newton–Euler).

## Results at a glance

All 12 cases are certified exactly (`certificates_summary.txt`).

| Robot | Joints | Base parameters, known dimensions | with rotor inertia and friction | Base parameters, symbolic dimensions |
|---|---|---|---|---|
| Universal Robots UR5e | 6 | 36 | 52 | 56 |
| ABB IRB 120 | 6 | 36 | 52 | 54 |
| Franka Emika Panda | 7 | 43 | 62 | 71 |

The known-dimension counts equal the closed-form value `1 + 7(n − 1)` for a revolute chain
with the first axis parallel to gravity. The explicit base-parameter lists are in
`<robot>/results/**/theta_base_list.txt` (and as LaTeX in `theta_base_latex.tex`).

See [RESULTS.md](RESULTS.md) for the full result tables and the index of result files.

## Requirements

Developed and tested with MATLAB R2026b and the Symbolic Math Toolbox. The Panda example
uses `loadrobot('frankaEmikaPanda')` (Robotics System Toolbox) only in `panda_robot.m`.

## Repository layout

```
UR5e/  IRB120/  Panda/
 ├── A_Full_Dynamics_<robot>_parallel_axis.m   symbolic Euler–Lagrange torques (parallel-axis inertia)
 ├── B_Creating_A_matrix_base.m                regrouping, symbolic link dimensions
 ├── B2_known_geometry_regrouping.m            regrouping, known link dimensions
 ├── B3/B4_actuation_*.m                       the same with rotor inertia and friction
 ├── regroup_base_parameters.m                 the five-stage framework (core function)
 ├── exact_verification.m                      polynomial-reduction certificate
 ├── torque_equations/                         generated torque equations (input to step B)
 └── results/                                  base-parameter lists and summaries
```

Panda only: `A2_regressor_columns_Panda.m` (monomial extraction by derivatives),
`panda_lagrange_torques.m`, `regressor_function.m` (block-wise regressor generation for large
models), `panda_rnea_numeric.m`.

## Quick start

```matlab
cd IRB120
B2_known_geometry_regrouping      % about 10–30 min; writes results/known_geometry/
```

The torque equations are included, so step A (symbolic derivation) can be skipped. To apply the
framework to **your own robot**, write the symbolic torque equations as MATLAB assignments
(`Torque_k = <expression>;`) in `torque_equations/Torque_1.txt … Torque_n.txt`, as in the examples
(see `load_torque_equations.m`), and adapt the kinematic symbols and settings at the top of the B scripts; the options of the core function are documented in `regroup_base_parameters.m`.

## Verification scripts

The saved regressors (`<robot>/results/**/base_regressor.mat`) let you rerun the exact checks without
repeating the regrouping:

| Script | What it proves | Output |
|---|---|---|
| `verification/certify_all.m` | exact dependency certificate `A_D = A_I*Gamma` for all 12 variants | `certificates_summary.txt` |
| `verification/reconstruction_check.m` | exact reconstruction `tau_sym = A*Theta` (six-joint robots; about 2–3 min per variant) | `reconstruction_check.txt` in each results folder |
| `verification/perturbation_test.m` | sampled residual versus exact verdict for a perturbed dependency coefficient (UR5e, IRB 120) | `perturbation_test.txt` |

```matlab
cd verification
V = reconstruction_check('UR5e', 'known_geometry');   % variants: results, known_geometry, actuation, actuation_known_geometry
```

The Panda reconstruction was checked numerically only (see the paper).
## Data and the identification study

The result summaries of the identification and validation study are included (see `RESULTS.md`). The Simscape datasets and the identification scripts that produced them are not part of
this lean release and can be provided on request.

## Citation

**If you use this software or the base-parameter results in your work, please cite the paper.**
GitHub's "Cite this repository" button (right-hand sidebar) provides the BibTeX and APA entries
from `CITATION.cff`. Until the paper is published, cite:

```bibtex
@unpublished{yusuf_base_parameters,
  author = {Yusuf, Muhammad Adel and Alyazidi, Nezar M. and Shafiullah, Md and Nasir, Ali},
  title  = {Automatic Symbolic--Numerical Derivation of Explicit Base Parameters from Robot Dynamic Equations},
  note   = {Submitted to Mechanism and Machine Theory. Code: https://github.com/Muhammad-Adel-Yusuf/base-parameter-regrouping},
  year   = {2026}
}
```

Archived at Zenodo: [10.5281/zenodo.23868741](https://doi.org/10.5281/zenodo.23868741). The DOI of the paper will be added when it is published.

## License and third-party material

Released under the [MIT License](LICENSE). The robot descriptions used in the
examples come from the ROS-Industrial (UR, ABB) and Franka Emika (`franka_description`,
Apache-2.0) projects; the Panda inertial data are from C. Gaz *et al.*, IEEE RA-L 4(4), 2019.

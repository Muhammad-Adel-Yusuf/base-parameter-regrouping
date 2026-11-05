# SYCIRO: SYmbolic Certified Identification of RObots

SYCIRO is a MATLAB tool for the dynamic modeling and identification of serial robots. From the
modified Denavit-Hartenberg (MDH) table of a robot it gives you

- the **symbolic dynamic model** (joint-torque equations, by Euler-Lagrange or Newton-Euler);
- the **base parameters**: the identifiable combinations of the inertial parameters, written explicitly
  in the physical symbols and checked by exact tests (dependency, reconstruction of the torque equations,
  minimality of the count);
- **least-squares identification** of the base parameters from your measured data or from a trajectory
  generated in the tool.

Revolute and prismatic joints are supported, link lengths can be known or left symbolic, and rotor inertia
and joint friction can be included. Each step can be used on its own; you do not have to run identification
to obtain the dynamic model or the base parameters.

## Requirements

- MATLAB with the **Symbolic Math Toolbox** (tested with R2026b)
- Optional: Parallel Computing Toolbox (faster runs)
- Only for the models in `simscape/`: Simulink, Simscape Multibody, and Robotics System Toolbox

## Quick start

1. Download the repository (**Code > Download ZIP**) and unzip it.
2. In MATLAB, open the unzipped folder and type

   ```matlab
   SYCIRO
   ```

3. In the **Robot** tab, choose a robot under *Load an example* (for instance the UR5e), click
   *Load example into the tables*, and press **RUN**.

## The interface

| Tab | What you do there |
|---|---|
| **Robot** | Enter the number of joints, the joint types, and the MDH table (numbers or expressions such as `pi/2`, `atan(3/4)`, `L1`, `L2+L3`). Give the link lengths (*Detect lengths* lists them; `NaN` = unknown). Choose known or symbolic lengths, rigid body or with rotor inertia and friction, and *What to compute*: **Base parameters and certificates**, or **Dynamic model only** (Euler-Lagrange or Newton-Euler). A robot typed here can be saved as `robots/mdh_<name>.m`. |
| **Results** | Certificate status of each variant, the list of base parameters, export as text, and the run log. |
| **Values** | Enter standard link data (mass, centre of mass, centroidal inertia, rotor inertia, friction) to obtain numerical base-parameter values. *Check against Newton-Euler* compares the model with an independent recursive Newton-Euler algorithm. |
| **Identify** | Least-squares identification for a variant with known lengths, either from your data (*Load data and identify*) or from an excitation trajectory generated in the tool, with optional noise and outliers. |

## Using your own data

The *Identify* tab reads

- a `.mat` file with four N-by-n arrays named `q`, `qd`, `qdd`, `tau` (one row per sample, one column per
  joint), or
- a `.csv` file with 4n columns `[q qd qdd tau]`.

Units are rad, rad/s, rad/s^2 and N m (m, m/s, m/s^2 and N for a prismatic joint). The joint coordinates
and torques must follow the MDH convention of the table you entered. If your excitation trajectory is a
MATLAB function, evaluate it at the sample times to obtain `q`, `qd`, `qdd`, add the measured torques, and
save the arrays:

```matlab
save('my_data.mat', 'q', 'qd', 'qdd', 'tau');
```

## Without the interface

```matlab
addpath mdh core robots
mdh = mdh_ur5e();                       % any table in robots/
out = base_parameters_from_mdh(mdh);    % all four variants
disp(struct2table(out.table))
```

The results are written to `results/<robot>/<variant>/`: `theta_base_list.txt` (the base parameters),
`regrouping_summary.txt` (ranks and residuals) and `exact_verification.txt` (the certificates). The four
variants are `symbolic` (lengths symbolic), `known` (lengths substituted), `actuation` and `actuation_known`
(with rotor inertia and friction). `examples/example_ur5e_base_parameters.m` goes step by step from the MDH
table to the base parameters, their numerical values, and their identification by least squares.

## Adding a robot

Create `robots/mdh_<name>.m`, or type the table in the interface and save it:

```matlab
function mdh = mdh_myrobot()
    syms L1 L2 real
    rows = [ 0     0   L1  0   0;       % [alpha_{i-1}  a_{i-1}  d_i  theta0_i  type], one row per joint
            -pi/2  0   0   0   0;       % type: 0 revolute, 1 prismatic
             0     L2  0   0   1];
    dims = struct('L1', 0.30, 'L2', 0.25);   % known values in metres
    mdh = mdh_table('MyRobot', rows, dims);
end
```

SYCIRO uses the modified DH convention of Khalil and Kleinfinger:
`{i-1}T{i} = Rot_x(alpha_{i-1}) Trans_x(a_{i-1}) Rot_z(theta_i) Trans_z(d_i)`, with `theta_i = q_i + theta0_i`
for a revolute joint and `d_i = q_i + d_i0` for a prismatic joint. Tables written in the classic (distal) DH
convention must be converted first.

## Included robots

| File | Robot |
|---|---|
| `mdh_ur5e` | Universal Robots UR5e |
| `mdh_irb120` | ABB IRB 120 |
| `mdh_panda` | Franka Emika Panda (7 joints) |
| `mdh_scara`, `mdh_scara_up` | SCARA (RRPR), prismatic axis drawn downward or upward |
| `mdh_stanford` | Stanford-type arm (RRPRRR) |
| `mdh_rpr_par`, `mdh_rpr_perp`, `mdh_rpr_twist` | RPR arms with parallel, perpendicular, and general-twist axes |

## Simscape models and simulated data

`simscape/<robot>/` contains, for the UR5e, IRB 120, and Panda:

- `simulink/<robot>_identification.slx`: Simscape Multibody model built from the robot description shipped
  with MATLAB (`simulink/generate_identification_data.m` rebuilds the model and regenerates the data);
- `data/simscape_data.mat`: identification trajectory (40 s at 250 Hz) and simulated joint torques;
- `data/validation_data.mat`: a separate 30 s trajectory for validation;
- `nominal_parameters_<robot>.m`: the reference inertial parameters in the MDH frames.

For the Panda, `simulink/view_panda_simulation.m` plays the motion in Mechanics Explorer.

The data files store the joint positions, velocities, accelerations, and torques in the MDH convention as
`t_data`, `td_data`, `tdd_data`, and `T_data` (N-by-n), with the sample times in `time_data`. To identify
from them in the *Identify* tab, save them under the names the tab expects:

```matlab
D = load('simscape/UR5e/data/simscape_data.mat');
q = D.t_data;  qd = D.td_data;  qdd = D.tdd_data;  tau = D.T_data;
save('ur5e_identification.mat', 'q', 'qd', 'qdd', 'tau');
```

## Scope

Serial chains with revolute and prismatic joints. The certificates concern the algebraic reduction of the
model defined by the entered table; they do not verify the table itself, the rigid-body assumptions, or the
real robot. The exact reconstruction check becomes slow above six joints and can be switched off.

## License and citation

MIT License (see `LICENSE`). The Panda link meshes come from the Franka Emika robot description
(Apache-2.0, see `simscape/Panda/simulink/meshes/LICENSE_franka_ros.txt`).

If you use SYCIRO, please cite:

> M. A. Yusuf, N. M. Alyazidi, M. Shafiullah, A. Nasir, *SYCIRO: From the Modified DH Table to Certified Base
> Parameters and Identification of Serial Robots*, submitted, 2026.

An earlier version of the software is archived on Zenodo: https://doi.org/10.5281/zenodo.23868741

# Key results

Summary of the results reported in the paper. All numbers are from simulated Simscape Multibody
data with known reference parameters; they verify the method and do not predict accuracy on hardware.

## 1. Base parameters (all 12 cases certified exactly)

| Robot | Joints | Dimensions | Without actuation terms | With rotor inertia and friction |
|---|---|---|---|---|
| UR5e | 6 | known | 36 | 52 |
| UR5e | 6 | symbolic | 56 | 72 |
| ABB IRB 120 | 6 | known | 36 | 52 |
| ABB IRB 120 | 6 | symbolic | 54 | 70 |
| Franka Emika Panda | 7 | known | 43 | 62 |
| Franka Emika Panda | 7 | symbolic | 71 | 90 |

- The known-dimension counts equal the closed-form value `1 + 7(n - 1)`.
- The structural rank was identical in all 60 numerical settings per case (3 state counts, 10 seeds, with
  and without column scaling); the smallest singular-value gap exceeded `1e13`.
- Each reduction is proved equivalent to the original torque equations by polynomial reduction
  modulo `sin^2 + cos^2 - 1` (see `certificates_summary.txt`).
- The identifiable subspace agrees with a numerical pivoted-QR reduction (largest principal angle
  about `1e-15` rad) and with an independent rank check on MATLAB's recursive Newton-Euler model.

## 2. Identification on the simulated trajectory (symbolic-dimension base sets: 56, 54 and 71 parameters)

| Robot | Estimator | Median relative error (%) | Torque RMSE (N m) |
|---|---|---|---|
| UR5e | OLS | 5.7e-8 | 6.5e-11 |
| UR5e | RLS | 3.1e-7 | 3.3e-10 |
| UR5e | Huber IRLS | 5.3e-8 | 6.5e-11 |
| IRB 120 | OLS | 2.1e-10 | 1.2e-14 |
| IRB 120 | RLS | 4.4e-6 | 1.0e-10 |
| IRB 120 | Huber IRLS | 2.1e-10 | 1.2e-14 |
| Panda | OLS | 1.5e-8 | 3.2e-11 |
| Panda | RLS | 4.1e-8 | 6.2e-11 |
| Panda | Huber IRLS | 1.8e-8 | 3.2e-11 |

## 3. Prediction on an independent validation trajectory (torque RMSE, N m)

| Robot | Torque RMS of the trajectory | OLS | RLS | Huber IRLS |
|---|---|---|---|---|
| UR5e | 11.87 | 1.3e-10 | 6.5e-10 | 1.3e-10 |
| IRB 120 | 6.86 | 1.2e-14 | 1.8e-10 | 1.2e-14 |
| Panda | 12.58 | 3.9e-11 | 9.0e-11 | 4.0e-11 |

## 4. Where to find the files

| What | File |
|---|---|
| Explicit base-parameter lists | `<robot>/results/**/theta_base_list.txt` and `theta_base_latex.tex` |
| Rank, residual and run time of each regrouping | `<robot>/results/**/regrouping_summary.txt` |
| Exact certificates | `certificates_summary.txt` |
| Sensitivity, numerical-QR baseline, independent check | `<robot>/results/benchmark_summary.txt` |
| Identification and validation summaries | `<robot>/results/identification_summary.txt`, `validation_summary.txt` |
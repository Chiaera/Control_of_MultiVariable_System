# LAB of Control of MultiVariable Systems

MATLAB exercises for a 3-DOF ROV (surge, sway, yaw), covering simulation, controllability, observability, and observer-based feedback control of the LTI model:

```
ẋ = Ax + Bu,
x = [x, y, θ, vx, vy, ω]ᵀ,
u = [Fx, Fy, Nz]ᵀ
```

## Requirements
- MATLAB with the **Symbolic Math Toolbox** (`syms`, `int`, `vpa` are used to compute the controllability Gramian in closed form).
- No external dependencies otherwise.

## Repository structure

```
Lab1/   Simulation of LTI systems
Lab2/   Controllability
Lab3/   Observability
Lab4/   Observer-based closed-loop control
```

---

### Lab1: Simulation of LTI systems
file: `Lab1.m` (+ `include/SimulateSysAnalytical.m`, `include/SimulateSysEuler.m`)

Simulates the ROV response to a constant input using both the closed-form analytical solution and forward-Euler integration, and compares the two.

Run: `Lab1.m` (no parameters to set).

---

### Lab2: Controllability
file: `Lab2.m`

Computes the reachability Gramian and the minimum-energy input steering `x0 = 0` to `xf = [-1, 2, 3, 0, 0, 0]` in `T = 20 s`, for three thruster configurations.

Set `scenario` at the top of the script:
| scenario | Case |
|---|---|
| 1 | All 4 thrusters enabled: fully controllable |
| 2 | Thruster 2 disabled: still fully controllable (redundancy) |
| 3 | Thrusters 2 & 3 disabled: **not** fully controllable (`rank(G) = 5`); prints the SVD (`U`, `D`, `V`) to identify the controllable subspace |

---

### Lab3: Observability
file: `Lab3.m`

Checks observability and runs a Luenberger observer **in open loop** (the nominal input from Lab2 is applied regardless of the estimate. Still **no feedback** here, that is Lab4). Plots the true vs. estimated state.

Set `scenario`:
| `scenario` | Case |
|---|---|
| `1` | Only position/yaw measured (`C` picks out `x, y, θ`) — dynamic Luenberger observer |
| `2` | Full state measured — trivial observer (`xHat = y`, no noise) |

---

### Lab4: Observer-based closed-loop control
file: `Lab4.m`

Same system, now with state feedback (`Kc`) correcting the nominal open-loop trajectory, optionally using an observer estimate instead of the true state, and optionally with noisy full-state measurements.

Set the three flags at the top of the script:
| Flags | Case |
|---|---|
| `enableClosedLoop=true, fullStateMeasured=flase, enableClosedLoopWithObserver=false` | **Q1:** full-state feedback |
| *(same flags as Q1, look at the estimation-error plot instead)* | **Q2:** Luenberger observer, position only, feedback still uses true state |
| `enableClosedLoop=true, enableClosedLoopWithObserver=true, fullStateMeasured=false` | **Q3:** feedback uses the Luenberger estimate |
| `enableClosedLoop=true, enableClosedLoopWithObserver=true, fullStateMeasured=true` | **Q4:** full state measured with sensor noise (`y_stdev=0.05`) |

**Notes:**
- The "Input" figure always plots the *nominal* feedforward input `uvalue`, not the actual applied `uref` (which includes the `-Kc·error` correction). It will look identical across Q1–Q4 by construction — this is expected, not a bug in the simulation itself, just something to be aware of when reading that particular plot.
- The Q3 observer gain (`wn=10`) is deliberately faster than the controller (`wn=1`); this causes a large but non-oscillatory transient overshoot in the estimation error ("peaking phenomenon" of high-gain observers) before it converges to ~0 within about 1–2 s.
- In the Q4 branch (`fullStateMeasured=true`), `xHat` is simply set equal to the noisy measurement at every step (no filtering). This means sensor noise passes straight through into the feedback gain `Kc` unattenuated, which is why Q4 looks noticeably jerkier than Q3. If a smoother Q4 result is desired, replace the trivial update with a real Luenberger correction (`Ko_loop` designed for `C = I`) instead of `xHat(:,k+1) = y_noise(:,k)`.

---

## Results
Each lab folder's `results/` subfolder contains the exported figures referenced in the write-up.

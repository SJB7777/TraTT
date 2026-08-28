# DynamicFTT — Program Design

Transient formation dynamics of acoustic Bragg satellites (KPS 2026 Fall abstract).

Design rule: **the solver produces `I(η,t)` and nothing else.** Every claim in the
paper is an analysis of that cube. No metric is added until a figure needs it.
This is the explicit fix for the previous repo becoming a black box.

---

## 0. Physics, fixed once

Two-beam Takagi–Taupin, one traveling acoustic modulation, stroboscopic
(quasi-static) X-ray response.

Justification for stroboscopic: the retardation group is `β = Ω·L/v_g ≈ Ω·L/c`.
At f = 1 GHz, L = 10 µm → β ≈ 2×10⁻⁴. Measured once in `tests/`, then dropped.
If a later result needs it, the channel-coupled form is a known extension.

### Fields

Normalized depth `ζ = z/L ∈ [0,1]`, slow time `t̃ = t/t_tr`, `t_tr = L/v_s`.

```
dD0/dζ = iπα · e^{+iφ(ζ;t̃)} · Dh
dDh/dζ = iπα · e^{-iφ(ζ;t̃)} · D0 + iη · Dh
```

### Drive

```
φ(ζ; t̃) = S · A(t̃ − ζ) · sin(q·ζ − w·t̃)
```

`A` = causal envelope, `A(s<0) = 0` gives the wavefront cutoff `z_max = v_s t`.

### Assumptions baked into "two free parameters"

`chi_0` is not assumed equal to anything — it is gauged away exactly. Substitute
`D_0 = e^{kz} d_0` and `D_h = e^{kz} d_h` with `k = i pi chi_0 / (lambda gamma)`,
the *same* factor for both beams; the `gamma k` term then cancels the `chi_0` term
in each equation. This works only because `chi_0` enters both equations with the
same coefficient, i.e. only in symmetric geometry (`gamma_0 = gamma_h`).

Be precise about complex `chi_0`: `Re chi_0` is a pure phase and does vanish from
`|D_h|^2`, but `Im chi_0` leaves a common `exp(-mu z)` that survives in the
absolute intensity. It cancels only in normalized observables — `p_n`, curve
shape, formation times — which is all we use. (This model sets `chi_0 = 0`
outright, so the point is moot here and matters only when comparing to data.)

`chi_0` does not vanish physically. It shifts the Bragg angle by
`|Re chi_0| / sin(2 theta_B)` — about 6.5 arcsec for Si(111) at 8 keV. That shift
is absorbed into the definition of `eta`, so `eta = 0` is the *refraction
corrected* Bragg angle, an offset to remember when matching real rocking data.

The real hidden assumption is in the coupling. It is properly
`C sqrt(chi_h chi_hbar)`, and writing it as `C|chi_h|` requires
`chi_hbar = chi_h*`, i.e. a **non-absorbing** crystal. That is exactly what
`alpha: float` in `core.py` encodes; with absorption `alpha` would be complex.

Full list:

1. two-beam
2. **symmetric geometry** `gamma_0 = gamma_h` — the `chi_0` gauge depends on it
3. plane wave, monochromatic
4. **no absorption**, `Im chi = 0`, so `alpha` is real
5. `chi_h chi_hbar = |chi_h|^2` (follows from 4; Si is centrosymmetric so the
   structure factor causes no extra trouble)
6. fixed polarization factor `C` (sigma: 1, pi: `cos 2 theta_B`), no mixing
7. the acoustic wave enters **only as a phase** — `chi_h` itself is unchanged
   (no elasto-optic change, no vibrational Debye–Waller reduction)

What 4 discards: uniform absorption (harmless — common factor) and the
**Borrmann effect** (not harmless — `Im chi_h` gives the two Bloch modes
different absorption, and mode-selective damping is itself a candidate extra
timescale). Deliberately excluded now, kept on the P8 mechanism list.

Note that strain is *not* discarded by 7: `d phi/d zeta = h . du/d zeta` is the
local lattice-spacing change, so it is already inside `phi`.

### Dimensionless groups — only two are free

| symbol | meaning | value |
|---|---|---|
| `α = L/Λ_ext` | dynamical diffraction strength | swept |
| `S = h·u₀` | acoustic modulation index | swept |
| `q = K·L` | carrier spatial phase across slab | `= 2πα_ac` |
| `w = Ω·t_tr` | acoustic phase per transit | `= q` (non-dispersive) |

`α_ac = L/Λ_ac` is the acoustic thickness parameter. **`q = w` is locked** by
`Ω = K v_s`; treating them as independent (previous repo) is a generalization we
do not need. Whether `α_ac` is tied to `α` or swept separately: see §6, decision D2.

### Analytic facts used as ground truth (not fitted)

1. Unmodulated symmetric Laue, η=0: `|Dh(1)|² = sin²(πα)`.
2. Lossless Laue: `|D0|² + |Dh|² = 1` at every ζ.
3. Phase grating expansion `e^{iφ} = Σ_n J_n(S) e^{in(qζ − w t̃)}` puts
   **channel n at `η_n = n·q`**. Satellite positions are known, never searched for.
4. Weak drive: `p_{±1}/p_0 → (J_1(S)/J_0(S))² → (S/2)²`.
5. Step envelope, t̃ ≥ 1: `φ` is already in its final form, so `I(η,t̃)` is exactly
   periodic in t̃ with period `2π/w`. **`t_f ≤ t_tr` is therefore structural.**
   Any run reporting `t_f > t_tr` under a step envelope is a bug, not a discovery.

Fact 5 reframes the central question of the abstract:

> not "does the transient outlive the wavefront?"
> but **"what fraction of the transit does formation actually take, and why?"**

### Standing prediction to test in Phase 4

Bragg geometry reflects only from the top ~`Λ_ext`. The front only has to cross
that depth for the diffracted state to be established:

```
t_f / t_tr  ≈  min(1, 1/α)      (Bragg)
t_f / t_tr  ≈  1                (Laue, whole thickness contributes)
```

Quantified in `explore/p0_equations.py`: for the **unmodulated static** problem,
the depth at which `|R|²` reaches 99 % of its final value is

| α | 0.5 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| ζ(99 %) | 0.983 | 0.866 | 0.480 | 0.244 | 0.127 | 0.060 |
| 1/α | 2.00 | 1.00 | 0.50 | 0.25 | 0.125 | 0.062 |

so `ζ_sat = min(1, 1/α)` holds to 2–4 % for α ≥ 4.

This is a static penetration depth, not a formation time — the real problem is
modulated and time dependent, so it may not carry over. That is exactly why it is
useful: a quantitative baseline fixed *before* the sweep, which P4 either
confirms or breaks. If it holds, the "additional dynamical timescale" of the
abstract is the extinction depth mapped into time, `t_ext = Λ_ext/v_s = t_tr/α`.

---

## 1. Module layout

Five files. Total target < 600 lines.

```
src/dftt/
  core.py      TT solver. frozen φ(ζ) → I(η). Laue (IVP) + Bragg (Riccati).
  drive.py     φ(ζ; t̃): traveling modulation, envelope, causality.
  cube.py      one parameter point → I(η,t̃) cube → .npz
  observe.py   I(η,t̃) → p_n(t), η_c,n(t), σ_n(t) → formation times
  sweep.py     (α,S) grid → runs.csv of t_f/t_tr
tests/
  test_core.py  the 5 analytic facts above, as asserts
```

Explicitly **not** ported from `../FTT`: `metrics.py`, `floquet.py`,
`experiment.py`, all of `scripts/`, all of `transient/`. They encode mechanism
hypotheses (quenching, Landau–Zener, asymmetry) that Phase 7 may or may not need.

`../AW_TT_EXP/spec_calc.py` stays where it is — it is the physical-units
feasibility calculator (Δθ, detector pixels, real t_tr) and is orthogonal to the
dimensionless solver. Called by hand when picking realistic (α,S).

### core.py

```python
def rocking(alpha, eta, phi, geometry="laue", n_z=512) -> I     # |Dh|² or |R|²
```

- **Laue**: initial-value problem, `D0(0)=1, Dh(0)=0`, RK4 forward, vectorized
  over the whole η grid at once.
- **Bragg**: two-point BVP. Solved as a Riccati IVP for `R = Dh/D0`, integrated
  backward from `ζ=1` with `R(1)=0`. One complex ODE, stable, no shooting.

`n_z` must resolve `q·ζ`: require `n_z ≥ 20·α_ac`. Asserted, not assumed.

### drive.py

Envelopes: `step` (default, exactly causal) and `smoothstep(rise)` for the finite
transducer rise time. `rise = T_r/t_tr` is the only physical knob that can push
formation past the transit under fact 5, so it is a deliberate, labeled control —
not a smoothing convenience.

### cube.py

One parameter point → one `.npz`:

```
eta  (Nη,) float64      # η grid, spans at least ±(M+0.5)·q
t    (Nt,) float64      # t̃ from 0 to t_max (default 3)
I    (Nt,Nη) float32    # the raw data
meta dict               # alpha, S, q, w, geometry, envelope, rise, n_z, M
```

Sampling: 32 t̃-samples per acoustic period `2π/w` → `Nt ≈ 32·(w/2π)·t_max`.
Nη = 512. float32 keeps a full sweep on disk in a few hundred MB; the derived
`p_n(t)` curves are tiny and always kept.

### observe.py

Given the cube:

- `p_n(t)` = ∫ I dη over `|η − n·q| < q/2` (window fixed analytically, §0 fact 3)
- `η_c,n(t)` = intensity centroid in that window
- `σ_n(t)` = second moment in that window

Formation time, **three independent definitions**, all reported:

| name | rule |
|---|---|
| `t_f^tol` | last t̃ at which `\|p_n − p̄_n^ss\| > ε·p̄_n^ss` (ε = 0.05) |
| `t_f^env` | decay time of the envelope of `p_n(t) − p̄_n^ss` |
| `t_f^int` | `∫\|p_n − p̄_n^ss\| dt / \|p_n(0) − p̄_n^ss\|` (no threshold at all) |

`p̄_n^ss` = mean over the last full acoustic period. Steady state here is a
**limit cycle, not a constant** — the drive keeps sliding after the front exits.
All three are computed for every run; if they disagree by more than ~2× the
result is flagged, not averaged.

### sweep.py

`(α, S)` grid → for each point run `cube` → `observe` → append a row to
`runs.csv`. Cubes for a labeled subset are kept; the rest keep only `p_n(t)`.
Multiprocessing over grid points.

---

## 2. Phases

Each phase ends with a concrete artifact. No phase starts before the previous
one's artifact exists.

| # | phase | artifact |
|---|---|---|
| P0 | `core.py` + `drive.py` + `tests/test_core.py` | 5 analytic asserts pass |
| P1 | `cube.py` | one `I(η,t̃)` movie, eyeballed |
| P2 | `observe.py` | `p_n(t)`, `η_c,n(t)`, `σ_n(t)` for that one cube |
| P3 | formation times | three `t_f` definitions agree on one cube |
| P4 | `sweep.py` over (α,S) | **`t_f/t_tr` map** — decision point |
| P5 | zoom on structure | which channel is slow? overshoot? ringing? |
| P6 | multiple timescales | `t_birth`, `t_peak`, `t_sat`, `t_relax` |
| P7 | scaling / collapse | does `t_f·α/t_tr` collapse? |
| P8 | mechanism | only now: extinction depth, Pendellösung, coupling |

P4 is the branch point named in the research plan. Its outcome decides P5–P8:

- map ≈ `min(1, 1/α)` → the extinction-depth story is the paper. Clean result.
- map has structure not explained by `1/α` → that structure is the paper.
- map ≈ 1 everywhere → formation is pure wavefront kinematics; the paper becomes
  the **formation sequence** (channel ordering `t_f,0 < t_f,1 < t_f,2`) instead
  of the formation time. Still a result, and still what the abstract promises.

All three outcomes are publishable. That is why the map comes before any mechanism.

---

## 3. What is deliberately not built

- No plotting library wrapper. matplotlib called directly in one `figs.py` per figure.
- No config system. Sweep parameters are literals at the top of `sweep.py`.
- No caching layer. `.npz` on disk is the cache.
- No abstract solver interface. Two geometries, two functions.
- No mechanism metrics (asymmetry, quenching, adiabaticity) until P8 asks for one.

---

## 4. Figure map (abstract → paper)

| fig | content | phase |
|---|---|---|
| 1 | `I(η,t̃)` movie / waterfall — transient exists | P1 |
| 2 | `p_0, p_±1, p_±2` vs t̃ — channel dynamics | P2 |
| 3 | characteristic times, several definitions | P3, P6 |
| 4 | `t_f/t_tr` over (α,S) | P4 |
| 5 | collapse onto a dimensionless time | P7 |
| 6 | mechanism | P8 |

---

## 5. Open decisions

**D1 — geometry.** Laue is an IVP (trivial); Bragg is a Riccati IVP (also easy)
and is where `t_f/t_tr < 1` is predicted to live. Recommend building both in P0;
the diff is ~30 lines and the Laue/Bragg contrast is Figure 4's control.

**D2 — is `α_ac` tied to `α`?** `α = L/Λ_ext` (X-ray) and `α_ac = L/Λ_ac`
(acoustic) are physically independent: `Λ_ext` is set by the reflection and
energy, `Λ_ac` by the acoustic frequency. Sweeping both is a 3D grid. Recommend
sweeping `(α, S)` at a few **fixed** `α_ac` values (e.g. 2, 5, 10) rather than a
full 3D grid — `α_ac` sets the satellite spacing `q`, and the plan's question is
about `α` and `S`.

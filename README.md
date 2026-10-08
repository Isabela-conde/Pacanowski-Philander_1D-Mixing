# Pacanowski–Philander 1-D Mixing

1-D vertical-mixing experiments with [OceanTurb.jl](https://github.com/glwagner/OceanTurb.jl), using
the **Pacanowski–Philander (1981)** Richardson-number scheme (and **KPP** for some runs). The repo
holds **two separate strands of work**, each in its own folder:

| folder | what it does | scripts |
|--------|--------------|---------|
| [`euc_shoaling/`](euc_shoaling) | how **shoaling the Equatorial Undercurrent** changes the surface flow through shear-driven mixing | `euc_pp.jl`, `strat_change_euc_pp.jl`, `pp_shear_driven.jl` |
| [`local_mixing/`](local_mixing) | eastern-Pacific warming: does surface warming **destabilise** the column and **warm the upwelled water**, and is it a **local** process? | `warming_mixing_test_pp.jl`, `upwelling_mixing_test_pp.jl`, `local_penetration_test_{pp,kpp}.jl`, `local_penetration_source_test_kpp.jl` |

Every figure lives in the same folder as the script that makes it, and each section below lists the
script → figure mapping. Run a script from inside its folder so the PNGs land next to it.

## The Pacanowski–Philander scheme

Each prognostic variable $\phi \in \{U, V, T, S\}$ obeys a 1-D vertical flux-divergence equation:

$$\frac{\partial \phi}{\partial t} = \frac{\partial}{\partial z}\left(K_\phi \, \frac{\partial \phi}{\partial z}\right) + R_\phi$$

where the only non-diffusive source $R_\phi$ is the Coriolis term on the horizontal velocities:

$$R_U = f\,V, \qquad R_V = -f\,U, \qquad R_T = R_S = 0$$

(At the equator $f = 0$, so the velocity equations reduce to pure vertical diffusion.) The
`local_mixing/` runs add further source terms to the $T$ equation — see that section.

### Richardson-number-dependent diffusivities

The eddy viscosity $K_U$ and diffusivity $K_T$ depend on the local gradient **Richardson number**:

$$Ri = \frac{N^2}{\left(\partial_z U\right)^2 + \left(\partial_z V\right)^2}, \qquad N^2 = \frac{\partial B}{\partial z}, \qquad B = g\left(\alpha T - \beta S\right)$$

$$K_U = \nu_0 + \frac{\nu_1}{\left(1 + c\,Ri\right)^{n}}, \qquad K_T = \kappa_0 + \frac{\kappa_1}{\left(1 + c\,Ri\right)^{n+1}}$$

with $K_V = K_U$ and $K_S = K_T$. As $Ri \to 0$ (strong shear / weak stratification) the diffusivities
approach their maxima $\nu_0 + \nu_1$, $\kappa_0 + \kappa_1$ — i.e. **shear instability switches mixing on**.
As $Ri \to \infty$ (strongly stratified) they relax to the background values $\nu_0$, $\kappa_0$.

| Symbol | Code | Value | Meaning |
|--------|------|-------|---------|
| $\nu_0$ | `Cν₀` | $10^{-4}$ | background viscosity |
| $\nu_1$ | `Cν₁` | $10^{-2}$ | max additional viscosity |
| $\kappa_0$ | `Cκ₀` | $10^{-5}$ | background diffusivity |
| $\kappa_1$ | `Cκ₁` | $10^{-2}$ | max additional diffusivity |
| $c$ | `Cc` | $5.0$ | Richardson coefficient |
| $n$ | `Cn` | $2.0$ | exponent |

---

## `euc_shoaling/` — EUC shoaling experiments

How the surface flow responds when the Equatorial Undercurrent is shoaled, through shear-driven
(Richardson-number) mixing.

### `euc_pp.jl`

- **Stratification** via a tanh thermocline temperature profile (S uniform), giving $N^2 = g\,\alpha\,\partial_z T$.
- **Initial velocity**: a westward surface jet (South Equatorial Current) over an eastward EUC at depth `h_core`:

$$U_0(z) = U_\text{sec}\,e^{z/h_\text{sec}} + U_\text{euc}\,\exp\!\left(-\frac{(-z - h_\text{core})^2}{2\,\sigma_\text{euc}^2}\right)$$

- **Comparison**: a *deep* EUC (`h_core = 130 m`) vs a *shoaled* EUC (`h_core = 80 m`), plus a sweep over core depth.
- **Forcing**: easterly surface wind stress applied as a flux boundary condition on $U$.

**Figures:** `euc_pp.png`, `euc_hovmoller.png`, `model_profiles.png`.

### `strat_change_euc_pp.jl` · `pp_shear_driven.jl`

`strat_change_euc_pp.jl` varies the stratification (thermocline strength) at fixed EUC and watches the
shear mixing respond → **figures** `strat_euc.png`, `strat_hovmoller.png`.
`pp_shear_driven.jl` is an `OceanTurb.jl` example modified to run with the Pacanowski–Philander scheme.

---

## `local_mixing/` — eastern-Pacific warming & upwelling tests

A sequence building toward one question: a talk claimed that under warming the eastern Pacific water
column becomes *less stable* and the *upwelled water warmer*, all by **local** (vertical) processes.
These runs test that, adding one ingredient at a time.

### Warming-sensitivity experiment (`warming_mixing_test_pp.jl`)

Tests a two-stage argument: *(1)* surface warming increases near-surface stratification and
**suppresses** mixing, then *(2)* continued warming weakens the thermocline stratification so $Ri$
drops back through $\tfrac14$ and mixing **returns**. The column keeps a fixed EUC+SEC shear and PP
mixing; a background $N^2$ and a surface heat-flux **feedback** keep it subcritical at $t=0$ and cap
the SST. Four runs isolate the two routes:

| Run | Forcing |
|-----|---------|
| `CTRL` | none |
| `SFC`  | surface heat flux (claim 1) |
| `THC`  | thermocline warming (claim 2) |
| `BOTH` | surface + thermocline |

**Key result:** surface heating alone stratifies the top and *traps* the heat — the thermocline $N^2$
is untouched and mixing never returns. Stage 2 happens **only** when the thermocline water itself is
warmed, which in reality is set remotely by the circulation, not by the local surface flux.

**Figures:** `warming_mixing_timeseries.png` (SST, near-surface & thermocline $N^2$, min $Ri$ in the
surface & thermocline bands, thermocline $\kappa$), `warming_mixing_hovmoller.png` ($Ri$ and $\kappa$
time–depth maps for `SFC` vs `THC`), `warming_mixing_profiles.png` ($T$ and $N^2$ profiles at start /
middle / end, all runs overlaid).

### Upwelling experiment (`upwelling_mixing_test_pp.jl`)

Adds the ingredient a pure mixing column cannot represent — **upwelling** — to test "warmer water was
upwelled". On top of the tamed setup it adds, all operator-split: a prescribed $w(z)>0$ (upwind
vertical advection of $T$), a deep **source-water** restoring toward $T_\text{profile}+\Delta T_\text{src}$,
the SST feedback, and a **convective adjustment** (see note). Runs: `NOUP`, `UPW`, `UPW+SFC` (local
surface warming), `UPW+SRC` (remote source warming), `UPW+BOTH`, plus a sweep over $w$.

**Key result:** SST is set by the upwelling rate (stronger $w$ → colder cold tongue). Warming SST has
two routes — *local* surface heating works without upwelling, but *remote* source warming reaches SST
**only when $w$ carries it up** ($\Delta\text{SST}\to 0$ as $w\to 0$): the "not local" point made quantitative.

**Figures:** `upwelling_mechanism.png`, `upwelling_sst_sensitivity.png`, `upwelling_anomaly_hovmoller.png`,
`upwelling_profiles.png`, `upwelling_penetration.png`, `upwelling_stability.png`, `upwelling_summary.png`.

> **Note (PP + static instability).** Because $K_T = \kappa_0 + \kappa_1/(1+c\,Ri)^{n+1}$ has an
> **odd** exponent, it goes *negative* once the column overturns ($N^2<0 \Rightarrow Ri<-1/c$),
> giving anti-diffusion and a numerical blow-up. Upwelling can create such inversions, so this script
> applies a convective adjustment each step (and tapers the source warming) to keep $N^2\ge0$.

### Local-penetration & source-water test (`local_penetration_test_pp.jl`, `local_penetration_test_kpp.jl`, `local_penetration_source_test_kpp.jl`)

These runs ask whether the destabilization *and* "warmer upwelled water" can be produced **locally** —
by surface heat penetrating down a vertical column — or whether the warm upwelled water has to be
**supplied remotely**. The `_pp` and `_kpp` scripts run the identical experiment with the two mixing
schemes; `_source_` adds a run that warms the source water directly.

#### Governing equations

Each column carries $U, V, T, S$ on a uniform grid ($N=128$, $H=300$ m); $S$ is fixed and the
equatorial ($f=0$) momentum equations are wind-forced and mixed as above. The **temperature** equation
collects every process explicitly:

$$\frac{\partial T}{\partial t} = \underbrace{\frac{\partial}{\partial z}\!\left(K_T\,\frac{\partial T}{\partial z}\right)}_{\text{(1) vertical mixing}} \;-\; \underbrace{w(z)\,\frac{\partial T}{\partial z}}_{\text{(2) upwelling}} \;+\; \underbrace{\frac{Q_\text{sw}}{\rho_0 c_p}\,\frac{e^{z/\lambda}}{\lambda}}_{\text{(3) penetrating shortwave}} \;+\; \underbrace{\frac{T^\star(z)-T}{\tau}\,\big[z\le z_\text{r}\big]}_{\text{(4) deep restoring / source}}$$

with the surface temperature flux (the **air–sea feedback**) as the top boundary condition:

$$Q_\theta = \frac{\lambda_\text{fb}}{\rho_0 c_p}\,\big(T_\text{sfc}-T_0\big) \qquad (\text{positive up: the ocean loses heat as it warms}).$$

| term | form | represents |
|------|------|-----------|
| (1) mixing | PP or KPP $K_T$ | shear-driven (+ boundary-layer) turbulence |
| (2) upwelling | $w(z)=w_0\,(1-e^{z/\delta_w})$, $\delta_w=15$ m, upwind | equatorial upwelling: 0 at the surface → $w_0$ below |
| (3) shortwave | Beer's law, e-folding $\lambda$; $\int dz = Q_\text{sw}$ | penetrating solar heating |
| (4) restoring | relax to $T^\star$ below $z_\text{r}$, $\tau=60$ d (the $[\cdot]$ is 1 only there) | the deep reservoir / supplied source water |
| BC | $Q_\theta=\lambda_\text{fb}(T_\text{sfc}-T_0)/\rho_0 c_p$ | air–sea damping that caps SST |

The restoring target is $T^\star(z) = T_\text{profile}(z) + \Delta T_\text{src}$ (ramped in over
$d_\text{ramp}$ so there is no step). With $\Delta T_\text{src}=0$ it is a cold **anchor** at depth
($z_\text{r}=-200$ m); with $\Delta T_\text{src}>0$ it **warms the upwelled source water**
($z_\text{r}=-100$ m in `_source_`, so the source depth is set directly) — the *remote* route. The
air–sea feedback is the key eastern-Pacific ingredient: it holds SST down while shortwave absorbed
*below* the mixed layer has no surface sink, so the subsurface can out-warm the damped surface.

#### Mixing scheme: PP vs KPP

- **PP** — interior Richardson-number mixing (as above); because $K_T$ goes negative for $N^2<0$, a
  **convective adjustment** is applied each step to keep $N^2\ge0$.
- **KPP** — the K-profile parameterization (Large–McWilliams–Doney 1994, as used in ACCESS-OM2): the
  surface buoyancy flux sets a boundary-layer depth with nonlocal transport, and static instability is
  removed internally (no adjustment needed). Terms (2)–(4) are passed through OceanTurb's native
  `Forcing`, and $Q_\theta$ is a function-of-model flux boundary condition so the feedback actually
  drives the boundary layer.

#### Experiments & figures

Two parameter regimes, each integrated 2 years with `:BackwardEuler`:

| regime | SW penetration $\lambda$ | upwelling $w_0$ | damping $\lambda_\text{fb}$ | heating $Q_\text{sw}$ |
|--------|--------------------------|-----------------|-----------------------------|------------------------|
| **favourable** | 50 m | 0.43 m/day | 50 W m⁻² K⁻¹ | 80 W m⁻² |
| **realistic** | 20 m | 1.5 m/day | 30 W m⁻² K⁻¹ | 80 W m⁻² |

`local_penetration_source_test_kpp.jl` adds the remote route, comparing four runs per regime: no heat;
SW only (local); source only (remote, $+1.5$ °C); SW + source.

| script | figures |
|--------|---------|
| `local_penetration_test_pp.jl` | `local_penetration_mechanism_pp.png`, `local_penetration_mechanism_pp_realistic.png`, `local_penetration_regime_pp.png` |
| `local_penetration_test_kpp.jl` | `local_penetration_mechanism_kpp.png`, `local_penetration_mechanism_kpp_realistic.png`, `local_penetration_regime_kpp.png` |
| `local_penetration_source_test_kpp.jl` | `local_source_kpp_favourable.png`, `local_source_kpp_realistic.png` |

**Key result.** With strong damping + deep penetration + weak upwelling, penetrating shortwave *can*
destabilize the **upper ~45 m** locally (KPP more than PP, because its boundary layer carries the heat
down), and the subsurface out-warms the damped surface. **But the source water (~120 m) barely warms** —
the anomaly is trapped in the boundary layer and decays with depth, and realistic penetration/upwelling
stabilizes instead. So the destabilization can be local, while **"warmer upwelled water" requires
warming the source directly (the remote route)** — the two cannot be produced together by vertical
processes alone.

## Acknowledgements

Developed with the assistance of **Claude** (Anthropic), which helped build diagnostics and Hovmöller plots, and document the model equations.

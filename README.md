# Pacanowski–Philander 1-D Mixing — EUC shoaling experiment

An [OceanTurb.jl](https://github.com/glwagner/OceanTurb.jl) experiment using the
**Pacanowski–Philander (1981)** Richardson-number-dependent vertical mixing scheme to test how
**shoaling the Equatorial Undercurrent (EUC)** changes the surface flow through shear-driven mixing, `euc_pp.jl` does this, whereas `pp_shear_driven.jl` is an `OceanTurb.jl` example modified to run with Pacanowski-Philander scheme.


## The Pacanowski–Philander scheme

Each prognostic variable $\phi \in \{U, V, T, S\}$ obeys a 1-D vertical flux-divergence equation:

$$\frac{\partial \phi}{\partial t} = \frac{\partial}{\partial z}\left(K_\phi \, \frac{\partial \phi}{\partial z}\right) + R_\phi$$

where the only non-diffusive source $R_\phi$ is the Coriolis term on the horizontal velocities:

$$R_U = f\,V, \qquad R_V = -f\,U, \qquad R_T = R_S = 0$$

(At the equator $f = 0$, so the velocity equations reduce to pure vertical diffusion.)

### Richardson-number-dependent diffusivities

The eddy viscosity $K_U$ and diffusivity $K_T$ depend on the local gradient **Richardson number**:

$$Ri = \frac{N^2}{\left(\partial_z U\right)^2 + \left(\partial_z V\right)^2}, \qquad N^2 = \frac{\partial B}{\partial z}, \qquad B = g\left(\alpha T - \beta S\right)$$

$$K_U = \nu_0 + \frac{\nu_1}{\left(1 + c\,Ri\right)^{n}}, \qquad K_T = \kappa_0 + \frac{\kappa_1}{\left(1 + c\,Ri\right)^{n+1}}$$

with $K_V = K_U$ and $K_S = K_T$. As $Ri \to 0$ (strong shear / weak stratification) the diffusivities
approach their maxima $\nu_0 + \nu_1$, $\kappa_0 + \kappa_1$ — i.e. **shear instability switches mixing on**.
As $Ri \to \infty$ (strongly stratified) they relax to the background values $\nu_0$, $\kappa_0$.

### Default parameters (PP81 / CV12)

| Symbol | Code | Value | Meaning |
|--------|------|-------|---------|
| $\nu_0$ | `Cν₀` | $10^{-4}$ | background viscosity |
| $\nu_1$ | `Cν₁` | $10^{-2}$ | max additional viscosity |
| $\kappa_0$ | `Cκ₀` | $10^{-5}$ | background diffusivity |
| $\kappa_1$ | `Cκ₁` | $10^{-2}$ | max additional diffusivity |
| $c$ | `Cc` | $5.0$ | Richardson coefficient |
| $n$ | `Cn` | $2.0$ | exponent |

## Experiment setup (`euc_pp.jl`)

- **Stratification** via a tanh thermocline temperature profile (S uniform), giving $N^2 = g\,\alpha\,\partial_z T$.
- **Initial velocity**: a westward surface jet (South Equatorial Current) over an eastward EUC at depth `h_core`:

$$U_0(z) = U_\text{sec}\,e^{z/h_\text{sec}} + U_\text{euc}\,\exp\!\left(-\frac{(-z - h_\text{core})^2}{2\,\sigma_\text{euc}^2}\right)$$

- **Comparison**: a *deep* EUC (`h_core = 130 m`) vs a *shoaled* EUC (`h_core = 80 m`), plus a sweep over core depth.
- **Forcing**: easterly surface wind stress applied as a flux boundary condition on $U$.


## Warming-sensitivity experiment (`warming_mixing_test_pp.jl`)

Tests a two-stage argument about the eastern tropical Pacific under warming: *(1)* surface
warming increases near-surface stratification and **suppresses** mixing, then *(2)* continued
warming weakens the thermocline stratification so $Ri$ drops back through $\tfrac14$ and mixing
**returns**. The column keeps a fixed EUC+SEC shear and PP mixing; a background $N^2$ and a surface
heat-flux **feedback** keep it subcritical at $t=0$ and cap the SST. Four runs isolate the two routes:

| Run | Forcing |
|-----|---------|
| `CTRL` | none |
| `SFC`  | surface heat flux (claim 1) |
| `THC`  | thermocline warming (claim 2) |
| `BOTH` | surface + thermocline |

**Key result:** surface heating alone stratifies the top and *traps* the heat — the thermocline
$N^2$ is untouched and mixing never returns. Stage 2 happens **only** when the thermocline water
itself is warmed, which in reality is set remotely by the circulation, not by the local surface flux.

Figures: `warming_mixing_timeseries.png` (SST, near-surface & thermocline $N^2$, min $Ri$ in the
surface & thermocline bands, thermocline $\kappa$), `warming_mixing_hovmoller.png` ($Ri$ and
$\kappa$ time–depth maps for `SFC` vs `THC`), and `warming_mixing_profiles.png` ($T$ and $N^2$
profiles at start / middle / end, all runs overlaid).

## Upwelling experiment (`upwelling_mixing_test_pp.jl`)

Adds the ingredient a pure mixing column cannot represent — **upwelling** — to test the last link,
"warmer water was upwelled". On top of the tamed setup it adds, all operator-split: a prescribed
$w(z)>0$ (upwind vertical advection of $T$), a deep **source-water** restoring toward
$T_\text{profile}+\Delta T_\text{src}$ (the remotely-supplied thermocline water), the SST feedback,
and a **convective adjustment** (see note below). Runs: `NOUP`, `UPW`, `UPW+SFC` (local surface
warming), `UPW+SRC` (remote source warming), `UPW+BOTH`, plus a sweep over $w$.

**Key result:** SST is set by the upwelling rate (stronger $w$ → colder cold tongue). Warming SST has
two routes — *local* surface heating works without upwelling, but *remote* source warming reaches SST
**only when $w$ carries it up** ($\Delta\text{SST}\to 0$ as $w\to 0$): the "not local" point made quantitative.

Figures: `upwelling_mechanism.png` (SST$(t)$, final $T(z)$, and `UPW+SRC` / `UPW+BOTH` anomaly
Hovmöllers), `upwelling_sst_sensitivity.png` (steady SST and $\Delta$SST vs $w$ for each warming
route), `upwelling_anomaly_hovmoller.png` ($T-T_\text{UPW}$ for `SFC` / `SRC` / `BOTH`), and
`upwelling_profiles.png` ($T$ and $N^2$ profiles at start / middle / end).

> **Note (PP + static instability).** Because $K_T = \kappa_0 + \kappa_1/(1+c\,Ri)^{n+1}$ has an
> **odd** exponent, it goes *negative* once the column overturns ($N^2<0 \Rightarrow Ri<-1/c$),
> giving anti-diffusion and a numerical blow-up. Upwelling can create such inversions, so this
> script applies a convective adjustment each step (and tapers the source warming) to keep $N^2\ge0$.

## Acknowledgements

Developed with the assistance of **Claude** (Anthropic), which helped build diagnostics and Hovmöller plots, and document the model equations.

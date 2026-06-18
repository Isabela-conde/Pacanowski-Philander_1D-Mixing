# Pacanowski–Philander 1-D Mixing — EUC shoaling experiment

A small set of [OceanTurb.jl](https://github.com/glwagner/OceanTurb.jl) experiments using the
**Pacanowski–Philander (1981)** Richardson-number-dependent vertical mixing scheme to test how
**shoaling the Equatorial Undercurrent (EUC)** changes the surface flow through shear-driven mixing.

> ⚠️ This is a 1-D column model with **no interior body force** (no zonal pressure gradient to hold
> the EUC up). The experiments are therefore **transient initial-value problems**: a jet structure is
> initialised and allowed to evolve under Ri-dependent mixing. The interesting physics lives in the
> early transient — see *Caveats* below.

## Contents

| File | Description |
|------|-------------|
| `pp_shear_driven.jl` | Minimal wind-driven shear-mixing example (reference for the OceanTurb PP API). |
| `euc_pp.jl` | EUC-shoaling sensitivity experiment: deep vs shoaled EUC, with Hovmöller diagnostics. |
| `euc_pp.png` | 3-panel summary: velocity profiles, surface u vs time, surface u vs EUC core depth. |
| `euc_hovmoller.png` | Time–depth Hovmöller of Ri and ν (eddy viscosity) for deep vs shoaled. |

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
- **Forcing**: optional easterly surface wind stress applied as a flux boundary condition on $U$.
- **Time stepping**: backward-Euler (unconditionally stable for the implicit diffusion).

## Running

```julia
# from the Julia REPL, in this directory
] add OceanTurb PyPlot     # one-time
include("euc_pp.jl")
```

Requires `OceanTurb`, `PyPlot` (and `ffmpeg` only if the optional U(z,t) animation block is enabled).
Outputs `euc_pp.png` and `euc_hovmoller.png`.

## Caveats — why long runs lose the signal

The velocity equation is **conservative vertical diffusion** with no momentum sink and no forcing that
sustains the jets. Two integral constraints follow:

- **Momentum is conserved**: $\frac{d}{dt}\int U\,dz =$ surface flux − bottom flux (zero in the unforced case).
- **Energy/shear is dissipated**: $\frac{d}{dt}\int \tfrac12 U^2\,dz = -\int \nu\,(\partial_z U)^2\,dz \le 0$.

So the eddy viscosity dissipates *shear*, and with nothing to maintain the EUC every initial condition
relaxes toward the **same well-mixed, depth-uniform profile** set by its (conserved) integrated momentum.
Because the deep and shoaled jets carry nearly the same integrated momentum, their profiles **converge at
long times**. This is a property of the experimental configuration, not a numerical instability — the
deep-vs-shoaled comparison is only meaningful during the early transient.

## Acknowledgements

Developed with the assistance of **Claude** (Anthropic), which helped adapt the OceanTurb PP API,
build the diagnostics and Hovmöller plots, and document the model equations and conservation properties.

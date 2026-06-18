# Pacanowski–Philander 1-D Mixing — EUC shoaling experiment

A small set of [OceanTurb.jl](https://github.com/glwagner/OceanTurb.jl) experiments using the
**Pacanowski–Philander (1981)** Richardson-number-dependent vertical mixing scheme to test how
**shoaling the Equatorial Undercurrent (EUC)** changes the surface flow through shear-driven mixing.


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



## Acknowledgements

Developed with the assistance of **Claude** (Anthropic), which helped build diagnostics and Hovmöller plots, and document the model equations.

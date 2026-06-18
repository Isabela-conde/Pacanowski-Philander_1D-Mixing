# EUC-shoaling sensitivity test using OceanTurb.jl's Pacanowski-Philander model.
#
# Idea (see notes): OceanTurb's PP model has NO interior body force (only
# Coriolis), so we can't hold the EUC up with a pressure gradient. Instead we
# run a TRANSIENT initial-value experiment: initialise a westward surface jet
# (SEC) over an eastward subsurface jet (EUC) at depth h_core, hold the
# stratification fixed via T(z), and let Ri-dependent mixing act. We compare how
# much the westward surface velocity erodes between a deep and a shoaled EUC.
#
# Stratification enters through temperature: N^2 = g*alpha*dT/dz (S uniform).
#
# Run:  julia euc_pp.jl        (needs: OceanTurb, PyPlot)
#       ] add OceanTurb PyPlot
#
# This file mirrors the API used in pp_shear_driven.jl (the working example):
#   - grid built with UniformGrid(N=, H=)
#   - Parameters imported from OceanTurb.PacanowskiPhilander
#   - initial conditions set by direct assignment (model.solution.U = f)
#   - boundary conditions built up-front and passed via bcs=
#   - time stepping with run_until!

using OceanTurb
using OceanTurb.PacanowskiPhilander: Parameters
using Printf
using PyPlot
using InteractiveUtils   # provides `varinfo`, which @use_pyplot_utils needs
                         # (auto-loaded in the REPL, but NOT in script mode)

@use_pyplot_utils

# ---------------------------------------------------------------------------
# Physical constants -- equator, so f = 0
# ---------------------------------------------------------------------------
constants = Constants(f = 0.0, α = 2.0e-4, β = 8.0e-4, ρ₀ = 1025.0, g = 9.81)

# ---------------------------------------------------------------------------
# Fixed stratification via a smooth thermocline temperature profile
#   T(z) warm at surface, cold at depth; N^2 = g*alpha*dT/dz peaks at the
#   thermocline. (z is NEGATIVE downward in OceanTurb.)
# ---------------------------------------------------------------------------
Tdeep, ΔT, h_th, δ_th = 14.0, 8.0, 30.0, 25.0
T₀(z) = Tdeep + 0.5ΔT * (1 + tanh((z + h_th) / δ_th))
S₀(z) = 35.0

# ---------------------------------------------------------------------------
# Two-jet initial velocity: westward SEC (surface) over eastward EUC at h_core
# ---------------------------------------------------------------------------
U_sec, h_sec = -0.30, 40.0      # westward surface jet [m/s], e-folding depth
U_euc, σ_euc = 0.80, 30.0       # EUC core speed [m/s], half-width
U₀(z, h_core) = U_sec * exp(z / h_sec) +
                U_euc * exp(-((-z - h_core)^2) / (2σ_euc^2))

# ---------------------------------------------------------------------------
# Build + run one experiment
# ---------------------------------------------------------------------------
function run_pp(h_core; N = 128, H = 300.0, mixing = true,
                Δt = 10minute, tfinal = 8day, wind_stress = 0.0)

    params = mixing ? Parameters() :                      # PP81/CV12 defaults
                      Parameters(Cν₁ = 0.0, Cκ₁ = 0.0)    # background only

    # Boundary conditions: optional easterly wind stress (flux on U at surface),
    # everything else insulating / no-flux. Built as FieldBoundaryConditions
    # (bottom, top) just like the working pp_shear_driven.jl example.
    Qᵘ = wind_stress / constants.ρ₀
    bcs = PacanowskiPhilander.BoundaryConditions(
        FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵘ)),  # U
        ZeroFluxBoundaryConditions(),                                                    # V
        ZeroFluxBoundaryConditions(),                                                    # T
        ZeroFluxBoundaryConditions(),                                                    # S
    )

    model = PacanowskiPhilander.Model(grid = UniformGrid(N = N, H = H),
                                      constants = constants,
                                      parameters = params,
                                      stepper = :BackwardEuler,
                                      bcs = bcs)

    # initial conditions (direct assignment of z-functions)
    model.solution.U = z -> U₀(z, h_core)
    model.solution.V = z -> 0.0
    model.solution.T = T₀
    model.solution.S = S₀

    # time stepping, sampling the surface velocity (top cell = index N)
    nsteps = Int(round(tfinal / Δt))
    usurf  = Float64[]; t = Float64[]
    push!(usurf, model.solution.U[N]); push!(t, 0.0)
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        push!(usurf, model.solution.U[N]); push!(t, n * Δt)
    end
    return model, t, usurf
end

# ---------------------------------------------------------------------------
# Experiments: deep vs shoaled EUC
# ---------------------------------------------------------------------------
h_deep, h_shoal = 130.0, 80.0
m_d, t, us_d = run_pp(h_deep)
m_s, _, us_s = run_pp(h_shoal)

@printf("surface u  deep EUC  (%.0f m): start %+.3f -> end %+.3f m/s\n",
        h_deep,  us_d[1], us_d[end])
@printf("surface u  shoaled   (%.0f m): start %+.3f -> end %+.3f m/s\n",
        h_shoal, us_s[1], us_s[end])

# core-depth sweep -> final surface velocity (the key result)
cores = 60.0:10.0:150.0
usurf_final = [run_pp(h)[3][end] for h in cores]

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
zc = nodes(m_d.solution.U)        # cell-centre depths

fig, ax = subplots(1, 3, figsize = (14, 5))

ax[1].plot([U₀(z, h_deep)  for z in zc], zc, "C0--", label = "deep init")
ax[1].plot([U₀(z, h_shoal) for z in zc], zc, "C3--", label = "shoaled init")
ax[1].plot([m_d.solution.U[i] for i in 1:length(zc)], zc, "C0", label = "deep final")
ax[1].plot([m_s.solution.U[i] for i in 1:length(zc)], zc, "C3", label = "shoaled final")
ax[1].axvline(0, color = "grey", lw = 0.7)
ax[1].set_xlabel("u [m/s]"); ax[1].set_ylabel("z [m]"); ax[1].legend(fontsize = 7)
ax[1].set_title("Velocity (dashed=init, solid=final)")

ax[2].plot(t ./ day, us_d, "C0", label = "deep EUC")
ax[2].plot(t ./ day, us_s, "C3", label = "shoaled EUC")
ax[2].set_xlabel("time [days]"); ax[2].set_ylabel("surface u [m/s]")
ax[2].legend(fontsize = 8); ax[2].set_title("Surface flow vs time")

ax[3].plot(collect(cores), usurf_final, "ko-")
ax[3].set_xlabel("EUC core depth [m]"); ax[3].set_ylabel("final surface u [m/s]")
ax[3].invert_xaxis(); ax[3].set_title("Surface flow vs EUC depth")

suptitle("OceanTurb PP: shoaling the EUC weakens the westward surface flow")
tight_layout()
savefig("euc_pp.png", dpi = 150)
println("saved euc_pp.png")

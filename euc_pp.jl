# EUC-shoaling sensitivity test using OceanTurb.jl's Pacanowski-Philander model.

using OceanTurb
using OceanTurb.PacanowskiPhilander: Parameters, KU, KT, local_richardson
using Printf
using PyPlot
using InteractiveUtils

@use_pyplot_utils

# Physical constants -- equator, so f = 0
constants = Constants(f = 0.0, α = 2.0e-4, β = 8.0e-4, ρ₀ = 1025.0, g = 9.81)

# stratification 
Tdeep, ΔT, h_th, δ_th = 14.0, 8.0, 30.0, 25.0
T₀(z) = Tdeep + 0.5ΔT * (1 + tanh((z + h_th) / δ_th))
S₀(z) = 35.0

# Two-jet initial velocity: westward SEC (surface) over eastward EUC at h_core
U_sec, h_sec = -0.30, 40.0      # westward surface jet [m/s], e-folding depth
U_euc, σ_euc = 0.80, 30.0       # EUC core speed [m/s], half-width
U₀(z, h_core) = U_sec * exp(z / h_sec) +
                U_euc * exp(-((-z - h_core)^2) / (2σ_euc^2))

#  thresholds
Ri_crit = 0.25
z_mld   = 70.0   

function pp_diagnostics(model)
    Ri = FaceField(model.grid)
    ν  = FaceField(model.grid)
    κ  = FaceField(model.grid)
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
        ν[i]  = KU(model, i)
        κ[i]  = KT(model, i)
    end
    return Ri, ν, κ
end

# Face profiles of Ri and ν at the current model state (plain Vectors), used to
# build the time-depth Hovmöller diagrams.
function snap_RiNu(model)
    Ri = FaceField(model.grid); ν = FaceField(model.grid)
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
        ν[i]  = KU(model, i)
    end
    return [Ri[i] for i in eachindex(Ri)], [ν[i] for i in eachindex(ν)]
end

# Minimum Ri inside the surface layer (finite faces only). A dip below Ri_crit
function min_Ri_upper(model; z_upper = z_mld)
    Ri = FaceField(model.grid)
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
    end
    zf   = nodes(Ri)
    vals = [Ri[i] for i in eachindex(Ri) if -z_upper <= zf[i] <= 0.0 && isfinite(Ri[i])]
    return isempty(vals) ? NaN : minimum(vals)
end

# Build and run one experiment
function run_pp(h_core; N = 128, H = 300.0, mixing = true,
                Δt = 10minute, tfinal = 6day, wind_stress = 0.05, 
                snapshots = false, save_every = 3)    # set wind stress

    params = mixing ? Parameters() :                      # PP81/CV12 defaults
                      Parameters(Cν₁ = 0.0, Cκ₁ = 0.0)    # background only

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

    # initial conditions 
    model.solution.U = z -> U₀(z, h_core)
    model.solution.V = z -> 0.0
    model.solution.T = T₀
    model.solution.S = S₀

    # time stepping; sample surface velocity (top cell = index N) and the
    # minimum Richardson number in the surface layer
    nsteps = Int(round(tfinal / Δt))
    usurf  = Float64[]; t = Float64[]; minRi = Float64[]
    tsnaps  = Float64[]
    Risnaps = Vector{Float64}[]; νsnaps = Vector{Float64}[]   # Ri(z,t), ν(z,t) Hovmöller
    zf = nodes(FaceField(model.grid))                  # face depths (shared by all snaps)

    function grab_snap()
        snapshots || return
        push!(tsnaps, model.clock.time)
        ri, nu = snap_RiNu(model)
        push!(Risnaps, ri); push!(νsnaps, nu)
    end

    push!(usurf, model.solution.U[N]); push!(t, 0.0)
    push!(minRi, min_Ri_upper(model)); grab_snap()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        push!(usurf, model.solution.U[N]); push!(t, n * Δt)
        push!(minRi, min_Ri_upper(model))
        snapshots && n % save_every == 0 && grab_snap()
    end
    return model, t, usurf, minRi, tsnaps, zf, Risnaps, νsnaps
end

# Experiments: deep vs shoaled EUC
h_deep, h_shoal = 130.0, 80.0
m_d, t, us_d, minRi_d, ts_d, zf, Ris_d, νs_d = run_pp(h_deep;  snapshots = true)
m_s, _, us_s, minRi_s, ts_s, _,  Ris_s, νs_s = run_pp(h_shoal; snapshots = true)

@printf("surface u  deep EUC  (%.0f m): start %+.3f -> end %+.3f m/s\n",
        h_deep,  us_d[1], us_d[end])
@printf("surface u  shoaled   (%.0f m): start %+.3f -> end %+.3f m/s\n",
        h_shoal, us_s[1], us_s[end])

# core-depth sweep -> final surface velocity vs EUC core depth
cores = 60.0:10.0:150.0
sweep       = [run_pp(h) for h in cores]
usurf_final = [s[3][end] for s in sweep]

# ------------------------------------------------------------------------------
# FIGURE 1 -- the original clean 3-panel summary
#   (1) velocity profiles   (2) surface u vs time   (3) core-depth sweep
# ------------------------------------------------------------------------------
zc = nodes(m_d.solution.U)        # cell-centre depths (U, T, ...)

fig, ax = subplots(1, 3, figsize = (14, 5))

# (1) velocity structure: dashed = initial, solid = final
ax[1].plot([U₀(z, h_deep)  for z in zc], zc, "C0--", label = "deep init")
ax[1].plot([U₀(z, h_shoal) for z in zc], zc, "C3--", label = "shoaled init")
ax[1].plot([m_d.solution.U[i] for i in 1:length(zc)], zc, "C0", label = "deep final")
ax[1].plot([m_s.solution.U[i] for i in 1:length(zc)], zc, "C3", label = "shoaled final")
ax[1].axvline(0, color = "grey", lw = 0.7)
ax[1].set_xlabel("u [m/s]"); ax[1].set_ylabel("z [m]"); ax[1].legend(fontsize = 7)
ax[1].set_title("Velocity")

# (2) surface velocity vs time
ax[2].plot(t ./ day, us_d, "C0", label = "deep EUC")
ax[2].plot(t ./ day, us_s, "C3", label = "shoaled EUC")
ax[2].set_xlabel("time [days]"); ax[2].set_ylabel("surface u [m/s]")
ax[2].legend(fontsize = 8); ax[2].set_title("Surface flow vs time")

# (3) core-depth sweep -> final surface velocity
ax[3].plot(collect(cores), usurf_final, "ko-")
ax[3].set_xlabel("EUC core depth [m]"); ax[3].set_ylabel("final surface u [m/s]")
ax[3].set_title("Surface flow vs EUC depth")

# suptitle("OceanTurb PP: shoaling the EUC weakens the westward surface flow")
tight_layout()
savefig("euc_pp.png", dpi = 150)
println("saved euc_pp.png")

# ------------------------------------------------------------------------------
# FIGURE 2 -- time-depth Hovmöller diagrams of Ri and ν, red contours of Ri_crit
# ------------------------------------------------------------------------------
# build [Nface, Ntime] matrices from the stored face-profile snapshots
hov(snaps) = reduce(hcat, snaps)          # each column is one time's face profile

RiM_d = hov(Ris_d);  RiM_s = hov(Ris_s)
νM_d  = hov(νs_d);   νM_s  = hov(νs_s)

LogNorm = PyPlot.matplotlib.colors.LogNorm
zmin    = -200.0                          # depth window to display

fig2, axmat2 = subplots(2, 2, figsize = (13, 8), sharex = true, sharey = true)
ax2 = vec(permutedims(axmat2))            # [Ri deep, Ri shoaled, ν deep, ν shoaled]

# --- Ri panels (clip to [0,1]; overlay the Ri_crit contour) -----------------
for (a, M, td, ttl) in ((ax2[1], RiM_d, ts_d, "Ri   deep EUC ($(Int(h_deep)) m)"),
                        (ax2[2], RiM_s, ts_s, "Ri   shoaled EUC ($(Int(h_shoal)) m)"))
    days = td ./ day
    pc = a.pcolormesh(days, zf, clamp.(M, 0, 1), cmap = "viridis",
                      vmin = 0, vmax = 1, shading = "auto")
    a.contour(days, zf, M, levels = [Ri_crit], colors = "r", linewidths = 1.3)
    a.set_ylim(zmin, 0); a.set_title(ttl); a.set_ylabel("z [m]")
    colorbar(pc, ax = a, label = "Ri")
end

# --- ν panels (log colour scale) -------------------------------------------
for (a, M, td, ttl) in ((ax2[3], νM_d, ts_d, "ν   deep EUC ($(Int(h_deep)) m)"),
                        (ax2[4], νM_s, ts_s, "ν   shoaled EUC ($(Int(h_shoal)) m)"))
    days = td ./ day
    pc = a.pcolormesh(days, zf, clamp.(M, 1e-6, 1e-1), cmap = "magma",
                      norm = LogNorm(vmin = 1e-5, vmax = 1e-2), shading = "auto")
    a.set_ylim(zmin, 0); a.set_title(ttl)
    a.set_xlabel("time [days]"); a.set_ylabel("z [m]")
    colorbar(pc, ax = a, label = "ν [m²/s]")
end

suptitle("Hovmöller diagrams of Ri and ν")
tight_layout()
savefig("euc_hovmoller.png", dpi = 150)
println("saved euc_hovmoller.png")


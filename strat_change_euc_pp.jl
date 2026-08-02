# EUC stratification-sensitivity test using OceanTurb.jl's Pacanowski-Philander model.
# Fix the EUC and vary the STRATIFICATION (thermocline strength), then watch how
# the EUC responds: stronger N² -> higher Ri -> less shear mixing -> EUC preserved;
# weaker N² -> lower Ri -> shear instability -> EUC eroded.

using OceanTurb
using OceanTurb.PacanowskiPhilander: Parameters, KU, KT, local_richardson
using Printf
using PyPlot
using InteractiveUtils

@use_pyplot_utils

# Physical constants -- equator, so f = 0
constants = Constants(f = 0.0, α = 2.0e-4, β = 8.0e-4, ρ₀ = 1025.0, g = 9.81)
# Thermocline shape is fixed; the STRATIFICATION is set per-run by ΔT (the
# temperature jump across the thermocline). Larger ΔT => stronger N².
Tdeep, h_th, δ_th = 14.0, 30.0, 25.0
Tprofile(z, ΔT) = Tdeep + 0.5ΔT * (1 + tanh((z + h_th) / δ_th))
S₀(z) = 35.0

# Peak buoyancy frequency at the thermocline for a given ΔT:  N² = g α ∂T/∂z,
# and ∂T/∂z is maximal at z = -h_th with value 0.5 ΔT / δ_th.
N2peak(ΔT) = constants.g * constants.α * 0.5ΔT / δ_th

# Fixed two-jet initial velocity: westward SEC (surface) over a fixed eastward EUC
h_core       = 100.0            # fixed EUC core depth [m]
U_sec, h_sec = -0.30, 40.0      # westward surface jet [m/s], e-folding depth
U_euc, σ_euc = 0.80, 30.0       # EUC core speed [m/s], half-width
U₀(z) = U_sec * exp(z / h_sec) +
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

# Face profiles of Ri and reduced shear squared, S²_red = (∂z U)² + (∂z V)² − 4N²,
# at the current model state (plain Vectors), used for the time-depth Hovmöllers.
# S²_red > 0  ⟺  Ri < 1/4  ⟺  shear-unstable.
function snap_RiSh(model)
    Ri = FaceField(model.grid); Sh = FaceField(model.grid)
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
        S² = ∂z(model.solution.U, i)^2 + ∂z(model.solution.V, i)^2
        N² = constants.g * (constants.α * ∂z(model.solution.T, i) -
                            constants.β * ∂z(model.solution.S, i))
        Sh[i] = S² - 4N²
    end
    return [Ri[i] for i in eachindex(Ri)], [Sh[i] for i in eachindex(Sh)]
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

# EUC core speed = peak eastward velocity in the column
u_core(model; N = model.grid.N) = maximum(model.solution.U[i] for i in 1:N)

# Build and run one experiment for a given stratification ΔT
function run_pp(ΔT; N = 128, H = 300.0, mixing = true,
                Δt = 10minute, tfinal = 6day, wind_stress = 0.01,
                snapshots = false, save_every = 3)

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

    # initial conditions: fixed velocity, stratification set by ΔT
    model.solution.U = z -> U₀(z)
    model.solution.V = z -> 0.0
    model.solution.T = z -> Tprofile(z, ΔT)
    model.solution.S = S₀

    # time stepping; sample EUC core speed, surface velocity (top cell = index N)
    # and the minimum Richardson number in the surface layer
    nsteps = Int(round(tfinal / Δt))
    usurf  = Float64[]; ucore = Float64[]; t = Float64[]; minRi = Float64[]
    tsnaps  = Float64[]
    Risnaps = Vector{Float64}[]; Shsnaps = Vector{Float64}[]   # Ri(z,t), S²_red(z,t) Hovmöller
    zf = nodes(FaceField(model.grid))                  # face depths (shared by all snaps)

    function grab_snap()
        snapshots || return
        push!(tsnaps, model.clock.time)
        ri, sh = snap_RiSh(model)
        push!(Risnaps, ri); push!(Shsnaps, sh)
    end

    push!(usurf, model.solution.U[N]); push!(ucore, u_core(model)); push!(t, 0.0)
    push!(minRi, min_Ri_upper(model)); grab_snap()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        push!(usurf, model.solution.U[N]); push!(ucore, u_core(model)); push!(t, n * Δt)
        push!(minRi, min_Ri_upper(model))
        snapshots && n % save_every == 0 && grab_snap()
    end
    return model, t, usurf, ucore, minRi, tsnaps, zf, Risnaps, Shsnaps
end

# Experiments: strong vs weak stratification
ΔT_strong, ΔT_weak = 10.0, 2.0
m_st, t, us_st, uc_st, mr_st, ts_st, zf, Ri_st, Sh_st = run_pp(ΔT_strong; snapshots = true)
m_wk, _, us_wk, uc_wk, mr_wk, ts_wk, _,  Ri_wk, Sh_wk = run_pp(ΔT_weak;   snapshots = true)

@printf("surface u  strong N² (ΔT=%.0f°C, N²=%.1e): start %+.3f -> end %+.3f m/s\n",
        ΔT_strong, N2peak(ΔT_strong), us_st[1], us_st[end])
@printf("surface u  weak   N² (ΔT=%.0f°C, N²=%.1e): start %+.3f -> end %+.3f m/s\n",
        ΔT_weak,   N2peak(ΔT_weak),   us_wk[1], us_wk[end])

# stratification sweep -> final surface velocity vs peak N²
strats      = 1.0:1.0:10.0                       # thermocline ΔT [°C]
sweep       = [run_pp(d) for d in strats]
usurf_final = [s[3][end] for s in sweep]
N2_sweep    = N2peak.(strats)

# ------------------------------------------------------------------------------
# FIGURE 1 -- 3-panel summary
#   (1) velocity profiles   (2) surface u vs time   (3) stratification sweep
# ------------------------------------------------------------------------------
zc = nodes(m_st.solution.U)        # cell-centre depths (U, T, ...)

fig, ax = subplots(1, 3, figsize = (14, 5))

# (1) velocity structure: initial (shared) dashed, final solid
ax[1].plot([U₀(z) for z in zc], zc, "k--", label = "initial")
ax[1].plot([m_st.solution.U[i] for i in 1:length(zc)], zc, "C0", label = "strong N² final")
ax[1].plot([m_wk.solution.U[i] for i in 1:length(zc)], zc, "C3", label = "weak N² final")
ax[1].axvline(0, color = "grey", lw = 0.7)
ax[1].set_xlabel("u [m/s]"); ax[1].set_ylabel("z [m]"); ax[1].legend(fontsize = 7)
ax[1].set_title("Velocity")

# (2) surface velocity vs time
ax[2].plot(t ./ day, us_st, "C0", label = "strong N² (ΔT=$(Int(ΔT_strong))°C)")
ax[2].plot(t ./ day, us_wk, "C3", label = "weak N² (ΔT=$(Int(ΔT_weak))°C)")
ax[2].set_xlabel("time [days]"); ax[2].set_ylabel("surface u [m/s]")
ax[2].legend(fontsize = 8); ax[2].set_title("Surface flow vs time")

# (3) stratification sweep -> final surface velocity
ax[3].plot(N2_sweep, usurf_final, "ko-")
ax[3].set_xlabel("peak N² [s⁻²]"); ax[3].set_ylabel("final surface u [m/s]")
ax[3].set_title("Surface flow vs stratification")

tight_layout()
savefig("strat_euc.png", dpi = 150)
println("saved strat_euc.png")

# ------------------------------------------------------------------------------
# FIGURE 2 -- time-depth Hovmöller diagrams of Ri and reduced shear squared
# ------------------------------------------------------------------------------
# build [Nface, Ntime] matrices from the stored face-profile snapshots
hov(snaps) = reduce(hcat, snaps)          # each column is one time's face profile

RiM_st = hov(Ri_st);  RiM_wk = hov(Ri_wk)
ShM_st = hov(Sh_st);  ShM_wk = hov(Sh_wk)

# symmetric colour scale for the (signed) reduced shear squared
shmax = maximum(abs.(filter(isfinite, vcat(vec(ShM_st), vec(ShM_wk)))))

zmin = -200.0                             # depth window to display

fig2, axmat2 = subplots(2, 2, figsize = (13, 8), sharex = true, sharey = true)
ax2 = vec(permutedims(axmat2))

# Ri
for (a, M, td, ttl) in ((ax2[1], RiM_st, ts_st, "Ri   strong N² (ΔT=$(Int(ΔT_strong))°C)"),
                        (ax2[2], RiM_wk, ts_wk, "Ri   weak N² (ΔT=$(Int(ΔT_weak))°C)"))
    days = td ./ day
    pc = a.pcolormesh(days, zf, clamp.(M, 0, 1), cmap = "viridis",
                      vmin = 0, vmax = 1, shading = "auto")
    a.contour(days, zf, M, levels = [Ri_crit], colors = "r", linewidths = 1.3)
    a.set_ylim(zmin, 0); a.set_title(ttl); a.set_ylabel("z [m]")
    colorbar(pc, ax = a, label = "Ri")
end

# reduced shear squared S²_red = S² − 4N²  (>0 ⟺ Ri < 1/4 ⟺ shear-unstable;
# black contour marks the marginal line S²_red = 0)
for (a, M, td, ttl) in ((ax2[3], ShM_st, ts_st, "S²_red   strong N² (ΔT=$(Int(ΔT_strong))°C)"),
                        (ax2[4], ShM_wk, ts_wk, "S²_red   weak N² (ΔT=$(Int(ΔT_weak))°C)"))
    days = td ./ day
    pc = a.pcolormesh(days, zf, M, cmap = "RdBu_r",
                      vmin = -shmax, vmax = shmax, shading = "auto")
    a.contour(days, zf, M, levels = [0.0], colors = "k", linewidths = 1.0)
    a.set_ylim(zmin, 0); a.set_title(ttl)
    a.set_xlabel("time [days]"); a.set_ylabel("z [m]")
    colorbar(pc, ax = a, label = "S²_red [s⁻²]")
end

suptitle("Hovmöller diagrams of Ri and reduced shear squared (S² − 4N²)")
tight_layout()
savefig("strat_hovmoller.png", dpi = 150)
println("saved strat_hovmoller.png")

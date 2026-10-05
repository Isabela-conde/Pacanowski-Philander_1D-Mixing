# warming_mixing_test_pp.jl
# -----------------------------------------------------------------------------
# Does the "warming -> restratify -> destratify -> mixing returns" story hold
# in a 1D column?  Same OceanTurb Pacanowski-Philander setup as
# strat_change_euc_pp.jl (equatorial column, f = 0, fixed EUC+SEC shear), but
# here we ADD time-evolving thermal forcing and watch the PP mixing respond.
#
# The argument is really two claims chained together:
#   (1) surface warming      -> stronger NEAR-SURFACE N²  -> Ri up   -> mixing suppressed
#   (2) thermocline warming  -> weaker  THERMOCLINE  N²  -> Ri < ¼  -> mixing returns
#                                                             (-> warm subsurface
#                                                                water reaches surface)
#
# The point of the test: claim (2) does NOT follow from claim (1) in a column.
# Surface heat stratifies the top and PP traps it there; the thermocline N² is
# untouched and the shear instability never comes back. You only get phase (2)
# if the thermocline WATER ITSELF warms. So we drive the two claims with their
# two different thermal forcings and see which produces what.
#
#   CTRL : no thermal forcing                     (baseline)
#   SFC  : surface heat flux only                 (tests claim 1; and whether 2 follows)
#   THC  : thermocline warming only               (tests claim 2)
#   BOTH : surface flux + thermocline warming      (the full "arc")
# -----------------------------------------------------------------------------

using OceanTurb
using OceanTurb.PacanowskiPhilander: Parameters, KU, KT, local_richardson
using Printf
using PyPlot

@use_pyplot_utils

# ---- physical constants (equator, f = 0) -----------------------------------
constants = Constants(f = 0.0, α = 2.0e-4, β = 8.0e-4, ρ₀ = 1025.0, g = 9.81)
cp = 3991.0   # seawater heat capacity [J kg⁻¹ K⁻¹]; only used to turn W m⁻² into a T flux

# ---- fixed thermocline shape; ΔT sets the INITIAL stratification ------------
Tdeep, h_th, δ_th = 14.0, 30.0, 25.0
Tprofile(z, ΔT) = Tdeep + 0.5ΔT * (1 + tanh((z + h_th) / δ_th))
S₀(z) = 35.0
N2peak(ΔT) = constants.g * constants.α * 0.5ΔT / δ_th   # peak N² at the thermocline

# ---- fixed two-jet shear: westward SEC over an eastward EUC -----------------
h_core       = 100.0            # EUC core depth [m]
U_sec, h_sec = -0.30, 40.0      # westward surface jet [m/s], e-folding depth
U_euc, σ_euc = 0.80, 30.0       # EUC core speed [m/s], half-width
U₀(z) = U_sec * exp(z / h_sec) + U_euc * exp(-((-z - h_core)^2) / (2σ_euc^2))

# ---- bands used to separate the two claims ----------------------------------
# Surface heating acts in the SURFACE band; the EUC shear / thermocline lives in
# the THERMOCLINE band. We diagnose min Ri and max κ in each band separately.
Z_SURF = (-40.0,   0.0)
Z_THC  = (-150.0, -40.0)
Ri_crit = 0.25
z_src   = -80.0                 # "upwelling source" depth: upper-thermocline water

# -----------------------------------------------------------------------------
# diagnostics
# -----------------------------------------------------------------------------
# min Ri and max κ (turbulent heat diffusivity) within a depth band [zlo, zhi]
function band_stats(model, zlo, zhi)
    Ri = FaceField(model.grid); κf = FaceField(model.grid)
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
        κf[i] = KT(model, i)
    end
    zf = nodes(Ri)
    minRi =  Inf; maxκ = -Inf
    for i in eachindex(Ri)
        if zlo <= zf[i] <= zhi
            isfinite(Ri[i]) && (minRi = min(minRi, Ri[i]))
            maxκ = max(maxκ, κf[i])
        end
    end
    return (isfinite(minRi) ? minRi : NaN), (isfinite(maxκ) ? maxκ : NaN)
end

# peak buoyancy frequency N² anywhere in the column [s⁻²]
function N2_peak_now(model)
    best = -Inf
    for i in eachindex(FaceField(model.grid))
        N² = constants.g * (constants.α * ∂z(model.solution.T, i) -
                            constants.β * ∂z(model.solution.S, i))
        best = max(best, N²)
    end
    return best
end

# temperature at the cell centre nearest a target depth
function T_at(model, ztarget)
    zc = nodes(model.solution.T)
    i  = argmin(abs.(zc .- ztarget))
    return model.solution.T[i]
end

# face profiles of Ri, reduced shear² (S² − 4N²), κ, and the cell T profile,
# for the time–depth Hovmöllers. S²_red > 0 ⟺ Ri < ¼ ⟺ shear-unstable.
function snap_profiles(model)
    Ri = FaceField(model.grid); Sh = FaceField(model.grid); κf = FaceField(model.grid)
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
        S² = ∂z(model.solution.U, i)^2 + ∂z(model.solution.V, i)^2
        N² = constants.g * (constants.α * ∂z(model.solution.T, i) -
                            constants.β * ∂z(model.solution.S, i))
        Sh[i] = S² - 4N²
        κf[i] = KT(model, i)
    end
    N = model.grid.N
    return ([Ri[i] for i in eachindex(Ri)],
            [Sh[i] for i in eachindex(Sh)],
            [κf[i] for i in eachindex(κf)],
            [model.solution.T[i] for i in 1:N])
end

# -----------------------------------------------------------------------------
# one experiment: fixed EUC+SEC shear, initial stratification ΔT, plus optional
# surface heat flux (Q_sfc, W m⁻², +ve = warming) and thermocline warming
# (r_thc, °C/day applied to the colder, deeper side of the thermocline).
# -----------------------------------------------------------------------------
function run_case(; label, ΔT = 8.0, Q_sfc = 0.0, r_thc = 0.0,
                    N = 128, H = 300.0, Δt = 10minute, tfinal = 120day,
                    wind_stress = 0.01, save_every = 144)

    # OceanTurb tracer flux is +ve upward (out of the ocean at the top), so a
    # surface heat INPUT of Q_sfc W m⁻² is a NEGATIVE top T flux. (Confirmed from
    # apply_top_bc!: rhs[N] -= flux/Δf, i.e. a -ve top flux warms the top cell.)
    Qᵀ = -Q_sfc / (constants.ρ₀ * cp)
    Qᵘ =  wind_stress / constants.ρ₀

    bcs = PacanowskiPhilander.BoundaryConditions(
        FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵘ)),  # U
        ZeroFluxBoundaryConditions(),                                                    # V
        FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵀ)),  # T
        ZeroFluxBoundaryConditions(),                                                    # S
    )

    model = PacanowskiPhilander.Model(grid = UniformGrid(N = N, H = H),
                                      constants = constants,
                                      parameters = Parameters(),     # PP81/CV12 defaults
                                      stepper = :BackwardEuler,
                                      bcs = bcs)

    model.solution.U = z -> U₀(z)
    model.solution.V = z -> 0.0
    model.solution.T = z -> Tprofile(z, ΔT)
    model.solution.S = S₀

    # thermocline-warming increment per step: deep-weighted shape (≈1 below the
    # thermocline, ≈0 at the surface) so warming the colder water SHRINKS the
    # vertical T contrast -> lowers N². Peak rate = r_thc °C/day.
    zc    = nodes(model.solution.T)
    gshp  = [0.5 * (1 - tanh((z + h_th) / δ_th)) for z in zc]
    dT_ps = r_thc * (Δt / day) .* gshp

    nsteps = Int(round(tfinal / Δt))
    t = Float64[]; sst = Float64[]; Tsrc = Float64[]
    N2mx = Float64[]
    minRi_s = Float64[]; maxκ_s = Float64[]
    minRi_t = Float64[]; maxκ_t = Float64[]
    tsnaps = Float64[]
    Risn = Vector{Float64}[]; Shsn = Vector{Float64}[]; κsn = Vector{Float64}[]; Tsn = Vector{Float64}[]
    zf = nodes(FaceField(model.grid))

    function record!()
        push!(t, model.clock.time)
        push!(sst, model.solution.T[N]); push!(Tsrc, T_at(model, z_src))
        push!(N2mx, N2_peak_now(model))
        rs, ks = band_stats(model, Z_SURF...); push!(minRi_s, rs); push!(maxκ_s, ks)
        rt, kt = band_stats(model, Z_THC... ); push!(minRi_t, rt); push!(maxκ_t, kt)
    end
    function snap!()
        push!(tsnaps, model.clock.time)
        ri, sh, kf, tt = snap_profiles(model)
        push!(Risn, ri); push!(Shsn, sh); push!(κsn, kf); push!(Tsn, tt)
    end

    record!(); snap!()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        @inbounds for i in 1:N            # apply thermocline warming (operator split)
            model.solution.T[i] += dT_ps[i]
        end
        record!()
        n % save_every == 0 && snap!()
    end

    @printf("[%-4s] ΔT₀=%.0f°C  Q_sfc=%+6.1f W/m²  r_thc=%.3f °C/day\n",
            label, ΔT, Q_sfc, r_thc)
    @printf("        SST      %6.3f -> %6.3f °C   (Δ = %+.2f)\n", sst[1], sst[end], sst[end]-sst[1])
    @printf("        T(%dm)   %6.3f -> %6.3f °C   (Δ = %+.2f)\n", Int(-z_src), Tsrc[1], Tsrc[end], Tsrc[end]-Tsrc[1])
    @printf("        peak N²  %.2e -> %.2e s⁻²\n", N2mx[1], N2mx[end])
    @printf("        surface  min Ri %.2f -> %.2f    max κ %.1e -> %.1e m²/s\n",
            minRi_s[1], minRi_s[end], maxκ_s[1], maxκ_s[end])
    @printf("        thermocl min Ri %.2f -> %.2f    max κ %.1e -> %.1e m²/s\n",
            minRi_t[1], minRi_t[end], maxκ_t[1], maxκ_t[end])

    return (; label, t, sst, Tsrc, N2mx, minRi_s, maxκ_s, minRi_t, maxκ_t,
              tsnaps, zf, zc, Risn, Shsn, κsn, Tsn)
end

# -----------------------------------------------------------------------------
# run the four experiments
# -----------------------------------------------------------------------------
Q0   = 120.0     # surface heating [W/m²]
Rthc = 0.06      # thermocline warming rate [°C/day]

ctrl = run_case(label = "CTRL", Q_sfc = 0.0,  r_thc = 0.0 )
sfc  = run_case(label = "SFC",  Q_sfc = Q0,   r_thc = 0.0 )
thc  = run_case(label = "THC",  Q_sfc = 0.0,  r_thc = Rthc)
both = run_case(label = "BOTH", Q_sfc = Q0,   r_thc = Rthc)
cases = (ctrl, sfc, thc, both)
cols  = Dict("CTRL"=>"k", "SFC"=>"C0", "THC"=>"C3", "BOTH"=>"C2")

# -----------------------------------------------------------------------------
# FIGURE 1 -- time series: does mixing switch off (claim 1) then on (claim 2)?
# -----------------------------------------------------------------------------
fig, ax = subplots(2, 2, figsize = (13, 9))

for c in cases
    d = c.t ./ day; col = cols[c.label]
    ax[1,1].plot(d, c.sst,     col, label = c.label)
    ax[1,2].plot(d, c.N2mx,    col, label = c.label)
    ax[2,1].plot(d, c.minRi_s, col, label = c.label)   # surface band -> claim 1
    ax[2,2].plot(d, c.minRi_t, col, label = c.label)   # thermocline band -> claim 2
end

ax[1,1].set_title("SST"); ax[1,1].set_ylabel("SST [°C]")
ax[1,2].set_title("peak N²  (stratification)"); ax[1,2].set_ylabel("N² [s⁻²]")
ax[2,1].set_title("min Ri, SURFACE band (0 to −40 m)  — claim 1")
ax[2,2].set_title("min Ri, THERMOCLINE band (−40 to −150 m)  — claim 2")
for a in (ax[2,1], ax[2,2])
    a.axhline(Ri_crit, color = "grey", ls = "--", lw = 1)
    a.text(0.02, Ri_crit + 0.01, "Ri = ¼ (shear-unstable below)", transform = a.get_yaxis_transform(), fontsize = 7, color = "grey")
    a.set_ylabel("min Ri"); a.set_ylim(0, 2)
end
for a in (ax[1,1], ax[1,2], ax[2,1], ax[2,2])
    a.set_xlabel("time [days]"); a.legend(fontsize = 8)
end
suptitle("Local warming → mixing response:  surface heating (claim 1) vs thermocline warming (claim 2)")
tight_layout()
savefig("warming_mixing_timeseries.png", dpi = 150)
println("saved warming_mixing_timeseries.png")

# -----------------------------------------------------------------------------
# FIGURE 2 -- Hovmöllers: SFC traps heat & keeps mixing off; THC restarts it
# -----------------------------------------------------------------------------
hov(snaps) = reduce(hcat, snaps)

fig2, axm = subplots(2, 3, figsize = (16, 8), sharex = true, sharey = true)
zmin = -200.0

function panel_T(a, c, ttl)
    pc = a.pcolormesh(c.tsnaps ./ day, c.zc, hov(c.Tsn), cmap = "inferno", shading = "auto")
    a.set_ylim(zmin, 0); a.set_title(ttl); a.set_ylabel("z [m]"); colorbar(pc, ax = a, label = "T [°C]")
end
function panel_Ri(a, c, ttl)
    M = hov(c.Risn)
    pc = a.pcolormesh(c.tsnaps ./ day, c.zf, clamp.(M, 0, 1), cmap = "viridis", vmin = 0, vmax = 1, shading = "auto")
    a.contour(c.tsnaps ./ day, c.zf, M, levels = [Ri_crit], colors = "r", linewidths = 1.3)
    a.set_ylim(zmin, 0); a.set_title(ttl); colorbar(pc, ax = a, label = "Ri")
end
function panel_κ(a, c, ttl)
    pc = a.pcolormesh(c.tsnaps ./ day, c.zf, log10.(max.(hov(c.κsn), 1e-6)), cmap = "magma", shading = "auto")
    a.set_ylim(zmin, 0); a.set_title(ttl); colorbar(pc, ax = a, label = "log₁₀ κ [m²/s]")
end

panel_T(axm[1,1], sfc, "SFC: T(z,t)")
panel_Ri(axm[1,2], sfc, "SFC: Ri  (red = ¼)")
panel_κ(axm[1,3], sfc, "SFC: mixing κ")
panel_T(axm[2,1], thc, "THC: T(z,t)")
panel_Ri(axm[2,2], thc, "THC: Ri  (red = ¼)")
panel_κ(axm[2,3], thc, "THC: mixing κ")
for a in (axm[2,1], axm[2,2], axm[2,3]); a.set_xlabel("time [days]"); end
suptitle("Surface heating (top) stratifies & traps → mixing stays off.  Thermocline warming (bottom) lowers Ri → mixing returns.")
tight_layout()
savefig("warming_mixing_hovmoller.png", dpi = 150)
println("saved warming_mixing_hovmoller.png")

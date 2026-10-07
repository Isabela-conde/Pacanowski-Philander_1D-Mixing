# local_penetration_test_pp.jl
# -----------------------------------------------------------------------------
# Can the speaker's story work LOCALLY?  Tests whether penetrating shortwave
# heating ALONE (no source-water warming, no horizontal advection) can warm the
# subsurface MORE than the surface and REDUCE stratification — i.e. make the
# column "less stable" by purely vertical, local processes.
#
# The enabling idea: the atmosphere damps SST hard (strong air-sea feedback), but
# shortwave deposited BELOW the mixed layer has no air-sea sink. If upwelling +
# mixing flush that subsurface heat out more slowly than it is deposited, the
# subsurface can out-warm the damped surface -> N² drops -> less stable. All local.
#
# Competition mapped here:  SW penetration depth λ (deposits heat deep)
#                     vs    upwelling w₀ (flushes it back up),  at strong damping.
#
# NOTE: same OceanTurb PP machinery as the other scripts (EUC+SEC shear, upwind
# upwelling of T, penetrating SW, SST feedback, convective adjustment). There is
# NO source-water warming: the deep restoring only anchors the abyss to a FIXED
# profile, so any subsurface warming must come from local SW penetration.
# -----------------------------------------------------------------------------

using OceanTurb
using OceanTurb.PacanowskiPhilander: Parameters, KU, KT, local_richardson
using Printf
using PyPlot

@use_pyplot_utils

constants = Constants(f = 0.0, α = 2.0e-4, β = 8.0e-4, ρ₀ = 1025.0, g = 9.81)
cp = 3991.0

Tdeep, h_th, δ_th = 14.0, 50.0, 30.0
N2_bg = 3.0e-5
Γ_bg  = N2_bg / (constants.g * constants.α)
Tprofile(z) = Tdeep + 0.5 * 10.0 * (1 + tanh((z + h_th) / δ_th)) + Γ_bg * z
S₀(z) = 35.0

h_core       = 110.0
U_sec, h_sec = -0.25, 50.0
U_euc, σ_euc = 0.60, 45.0
U₀(z) = U_sec * exp(z / h_sec) + U_euc * exp(-((-z - h_core)^2) / (2σ_euc^2))

δw = 15.0
w_prof(z, w₀) = w₀ * (1 - exp(z / δw))

z_restore = -200.0        # anchor only the abyss (deep, fixed)
τ_src     = 60day
z_sub     = -50.0         # "subsurface" diagnostic depth (thermocline)

function convect!(T, N; maxsweeps = 200)
    for _ in 1:maxsweeps
        stable = true
        @inbounds for i in 2:N
            if T[i] < T[i-1]
                m = 0.5 * (T[i] + T[i-1]); T[i] = m; T[i-1] = m; stable = false
            end
        end
        stable && break
    end
end

function T_at(model, ztarget)
    zc = nodes(model.solution.T); i = argmin(abs.(zc .- ztarget))
    return model.solution.T[i]
end

# mean N² over a depth band from the current model state
function meanN2(model, zlo, zhi)
    T = model.solution.T; zc = nodes(T); n = length(zc)
    s = 0.0; c = 0
    @inbounds for i in 2:n-1
        if zlo <= zc[i] <= zhi
            s += constants.g * constants.α * (T[i+1] - T[i-1]) / (zc[i+1] - zc[i-1]); c += 1
        end
    end
    return c == 0 ? NaN : s / c
end

# N² profile from a cell-T snapshot (for the final-profile panels)
function n2_profile(Tcol, zc)
    n = length(zc); zmid = zc[2:n-1]
    N2 = [constants.g * constants.α * (Tcol[i+1] - Tcol[i-1]) / (zc[i+1] - zc[i-1]) for i in 2:n-1]
    return zmid, N2
end

# -----------------------------------------------------------------------------
# one run: penetrating SW heating Q_sfc (e-folding λ_pen), SST damping λ_fb,
# upwelling w₀.  No source warming.
# -----------------------------------------------------------------------------
function run_case(; Q_sfc = 0.0, λ_pen = 30.0, λ_fb = 50.0, w₀ = 1.0e-5,
                    N = 128, H = 300.0, Δt = 20minute, tfinal = 365day,
                    wind_stress = 0.01, rec_every = 72)

    Qᵘ = wind_stress / constants.ρ₀
    bcs = PacanowskiPhilander.BoundaryConditions(
        FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵘ)),  # U
        ZeroFluxBoundaryConditions(),                                                    # V
        ZeroFluxBoundaryConditions(),                 # T: no surface flux (heating is SW)
        ZeroFluxBoundaryConditions(),                                                    # S
    )
    model = PacanowskiPhilander.Model(grid = UniformGrid(N = N, H = H), constants = constants,
                                      parameters = Parameters(), stepper = :BackwardEuler, bcs = bcs)
    model.solution.U = z -> U₀(z); model.solution.V = z -> 0.0
    model.solution.T = z -> Tprofile(z); model.solution.S = S₀

    Δz   = H / N
    zc   = nodes(model.solution.T)
    SST₀ = model.solution.T[N]
    wz   = [w_prof(z, w₀) for z in zc]
    rad  = Q_sfc > 0 ? [(Q_sfc / (constants.ρ₀ * cp)) * exp(z / λ_pen) / λ_pen * Δt for z in zc] :
                       zeros(length(zc))
    Ttgt = [Tprofile(z) for z in zc]                 # fixed deep anchor (no warming)

    nsteps = Int(round(tfinal / Δt))
    t = Float64[]; sst = Float64[]; Tsub = Float64[]; n2i = Float64[]
    rec!() = (push!(t, model.clock.time); push!(sst, model.solution.T[N]);
              push!(Tsub, T_at(model, z_sub)); push!(n2i, meanN2(model, -100.0, 0.0)))

    rec!()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        @inbounds for i in N:-1:2                     # upwelling advection of T
            model.solution.T[i] += -wz[i] * (model.solution.T[i] - model.solution.T[i-1]) / Δz * Δt
        end
        if Q_sfc > 0                                  # penetrating shortwave
            @inbounds for i in 1:N; model.solution.T[i] += rad[i]; end
        end
        @inbounds for i in 1:N                        # deep anchor (abyss only)
            zc[i] <= z_restore && (model.solution.T[i] += (Δt / τ_src) * (Ttgt[i] - model.solution.T[i]))
        end
        model.solution.T[N] -= λ_fb * (model.solution.T[N] - SST₀) * Δt / (constants.ρ₀ * cp * Δz)  # SST damping
        convect!(model.solution.T, N)
        n % rec_every == 0 && rec!()
    end
    rec!()                                            # capture final state
    Tfinal = [model.solution.T[i] for i in 1:N]
    return (; t, sst, Tsub, n2i, zc, Tfinal, n2ss = n2i[end])
end

# -----------------------------------------------------------------------------
# A -- illustrative run: strong damping, deep SW, weak upwelling, 2 years
# -----------------------------------------------------------------------------
w0_A, λfb_A, λpen_A, Q_A = 0.5e-5, 50.0, 50.0, 80.0
baseA = run_case(Q_sfc = 0.0, w₀ = w0_A, λ_fb = λfb_A, Δt = 15minute, tfinal = 730day, rec_every = 96)
heatA = run_case(Q_sfc = Q_A, λ_pen = λpen_A, w₀ = w0_A, λ_fb = λfb_A, Δt = 15minute, tfinal = 730day, rec_every = 96)

@printf("illustrative run (λ_pen=%.0fm, w₀=%.1f m/day, λ_fb=%.0f):\n", λpen_A, w0_A*86400, λfb_A)
@printf("  final ΔSST = %+.2f °C,  ΔT(%dm) = %+.2f °C\n",
        heatA.sst[end]-baseA.sst[end], Int(-z_sub), heatA.Tsub[end]-baseA.Tsub[end])
@printf("  ⟨N²⟩ 0–100 m: no-heat %.2e -> heated %.2e  (Δ = %+.2e)\n",
        baseA.n2ss, heatA.n2ss, heatA.n2ss - baseA.n2ss)

figA, ax = subplots(2, 2, figsize = (13, 9))

ax[1,1].plot(heatA.t ./ day, heatA.sst  .- baseA.sst,  "C3", label = "surface (SST)")
ax[1,1].plot(heatA.t ./ day, heatA.Tsub .- baseA.Tsub, "C0", label = @sprintf("subsurface (%d m)", Int(-z_sub)))
ax[1,1].set_xlabel("time [days]"); ax[1,1].set_ylabel("ΔT [°C]")
ax[1,1].set_title("Warming: surface vs subsurface"); ax[1,1].legend(fontsize = 8)

ax[1,2].plot(baseA.t ./ day, baseA.n2i, "k",  label = "no heating")
ax[1,2].plot(heatA.t ./ day, heatA.n2i, "C3", label = "penetrating SW")
ax[1,2].set_xlabel("time [days]"); ax[1,2].set_ylabel("⟨N²⟩ 0–100 m [s⁻²]")
ax[1,2].set_title("Upper-ocean stratification"); ax[1,2].legend(fontsize = 8)

ax[2,1].plot(heatA.Tfinal .- baseA.Tfinal, heatA.zc, "C3")
ax[2,1].axvline(0, color = "grey", lw = 0.6); ax[2,1].set_ylim(-150, 0)
ax[2,1].set_xlabel("ΔT [°C]"); ax[2,1].set_ylabel("z [m]"); ax[2,1].set_title("Final ΔT (heated − no heat)")

zmid, nh = n2_profile(heatA.Tfinal, heatA.zc)
_,    nb = n2_profile(baseA.Tfinal, baseA.zc)
ax[2,2].plot(nh .- nb, zmid, "C3")
ax[2,2].axvline(0, color = "grey", lw = 0.6); ax[2,2].set_ylim(-150, 0)
ax[2,2].set_xlabel("ΔN² [s⁻²]"); ax[2,2].set_title("Final ΔN² (heated − no heat)")

suptitle("Local penetrating shortwave under strong SST damping")
tight_layout()
savefig("local_penetration_mechanism.png", dpi = 150)
println("saved local_penetration_mechanism.png")

# -----------------------------------------------------------------------------
# B -- regime map: does the column destratify as a function of (λ_pen, w₀)?
#      ΔN² = ⟨N²⟩(heated) − ⟨N²⟩(no heat).  Negative = destratified = local works.
# -----------------------------------------------------------------------------
w0s   = collect(0.0:0.5:3.0) .* 1e-5      # upwelling [m/s]
λpens = collect(10.0:10.0:60.0)           # SW e-folding [m]
ΔN2   = fill(NaN, length(λpens), length(w0s))

for (j, w) in enumerate(w0s)
    base = run_case(Q_sfc = 0.0, w₀ = w, λ_fb = 50.0, Δt = 30minute, tfinal = 365day)
    for (i, λ) in enumerate(λpens)
        h = run_case(Q_sfc = 80.0, λ_pen = λ, w₀ = w, λ_fb = 50.0, Δt = 30minute, tfinal = 365day)
        ΔN2[i,j] = h.n2ss - base.n2ss
    end
    @printf("w₀=%.1f m/day done\n", w * 86400)
end

figB, bx = subplots(1, 1, figsize = (7.5, 6))
amax = maximum(abs.(filter(isfinite, vec(ΔN2))))
pc = bx.pcolormesh(w0s .* 86400, λpens, ΔN2, cmap = "RdBu_r", vmin = -amax, vmax = amax, shading = "auto")
bx.contour(w0s .* 86400, λpens, ΔN2, levels = [0.0], colors = "k", linewidths = 1.2)
bx.set_xlabel("upwelling w₀ [m/day]"); bx.set_ylabel("SW penetration depth λ [m]")
bx.set_title("Δ⟨N²⟩ 0–100 m  (heated − no heat)")
colorbar(pc, ax = bx, label = "Δ⟨N²⟩ [s⁻²]")
tight_layout()
savefig("local_penetration_regime.png", dpi = 150)
println("saved local_penetration_regime.png")

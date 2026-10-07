# local_penetration_test_kpp.jl
# -----------------------------------------------------------------------------
# KPP version of local_penetration_test_pp.jl -- same experiment, swapped mixing
# scheme (OceanTurb.KPP, Large-McWilliams-Doney 1994 K-profile parameterization).
#
# Same question: can LOCAL penetrating shortwave (no source-water warming, no
# horizontal advection) warm the subsurface more than the surface and reduce
# stratification?  Run this against the PP version to see whether the answer is
# scheme-robust.
#
# KPP differences wired in here:
#   * the surface heat flux drives the KPP boundary layer -> the air-sea feedback
#     is given as a FUNCTION flux BC (Qθ depends on SST), not an operator-split.
#   * penetrating shortwave, upwelling advection and the deep anchor are passed
#     through KPP's native Forcing interface (tendencies, per second).
#   * NO convective adjustment: KPP removes static instability internally.
# -----------------------------------------------------------------------------

using OceanTurb
const KPP = OceanTurb.KPP        # qualify to avoid a Parameters/Model name clash
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

z_restore = -200.0
τ_src     = 60day
z_sub     = -50.0

function T_at(model, ztarget)
    zc = nodes(model.solution.T); i = argmin(abs.(zc .- ztarget))
    return model.solution.T[i]
end

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

function n2_profile(Tcol, zc)
    n = length(zc); zmid = zc[2:n-1]
    N2 = [constants.g * constants.α * (Tcol[i+1] - Tcol[i-1]) / (zc[i+1] - zc[i-1]) for i in 2:n-1]
    return zmid, N2
end

# -----------------------------------------------------------------------------
# one KPP run
# -----------------------------------------------------------------------------
function run_case(; Q_sfc = 0.0, λ_pen = 30.0, λ_fb = 50.0, w₀ = 1.0e-5,
                    N = 128, H = 300.0, Δt = 20minute, tfinal = 365day,
                    wind_stress = 0.01, rec_every = 72)

    grid = UniformGrid(N, H)
    zc   = grid.zc
    Δz   = H / N
    SST₀ = Tprofile(zc[N])
    wz   = [w_prof(z, w₀) for z in zc]
    rad  = Q_sfc > 0 ? [(Q_sfc / (constants.ρ₀ * cp)) * exp(z / λ_pen) / λ_pen for z in zc] :
                       zeros(length(zc))          # SW heating RATE [°C/s]
    Ttgt = [Tprofile(z) for z in zc]

    # T forcing (tendency, per second): penetrating SW + upwind upwelling of T +
    # deep anchor restoring. (KPP integrates this explicitly via R_T.)
    function forcing_T(m, i)
        f = rad[i]
        i > 1        && (f += -wz[i] * (m.solution.T[i] - m.solution.T[i-1]) / Δz)
        zc[i] <= z_restore && (f += (Ttgt[i] - m.solution.T[i]) / τ_src)
        return f
    end

    # air-sea feedback: surface T flux (+ve up = heat leaving) that drives the KPP
    # boundary layer and damps SST toward SST₀.
    Qθ(m) = λ_fb * (m.solution.T[N] - SST₀) / (constants.ρ₀ * cp)
    Qᵘ    = wind_stress / constants.ρ₀

    bcs = KPP.ModelBoundaryConditions(
        U = FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵘ)),
        T = FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qθ)),
    )

    model = KPP.Model(grid = grid, constants = constants, parameters = KPP.Parameters(),
                      stepper = :BackwardEuler, bcs = bcs, forcing = KPP.Forcing(T = forcing_T))

    model.solution.U = z -> U₀(z); model.solution.V = z -> 0.0
    model.solution.T = z -> Tprofile(z); model.solution.S = S₀

    nsteps = Int(round(tfinal / Δt))
    t = Float64[]; sst = Float64[]; Tsub = Float64[]; n2i = Float64[]
    rec!() = (push!(t, model.clock.time); push!(sst, model.solution.T[N]);
              push!(Tsub, T_at(model, z_sub)); push!(n2i, meanN2(model, -100.0, 0.0)))

    rec!()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        n % rec_every == 0 && rec!()
    end
    rec!()
    Tfinal = [model.solution.T[i] for i in 1:N]
    return (; t, sst, Tsub, n2i, zc, Tfinal, n2ss = n2i[end])
end

# -----------------------------------------------------------------------------
# A -- illustrative run: strong damping, deep SW, weak upwelling, 2 years
# -----------------------------------------------------------------------------
w0_A, λfb_A, λpen_A, Q_A = 0.5e-5, 50.0, 50.0, 80.0
baseA = run_case(Q_sfc = 0.0, w₀ = w0_A, λ_fb = λfb_A, Δt = 15minute, tfinal = 730day, rec_every = 96)
heatA = run_case(Q_sfc = Q_A, λ_pen = λpen_A, w₀ = w0_A, λ_fb = λfb_A, Δt = 15minute, tfinal = 730day, rec_every = 96)

@printf("KPP illustrative run (λ_pen=%.0fm, w₀=%.1f m/day, λ_fb=%.0f):\n", λpen_A, w0_A*86400, λfb_A)
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

suptitle("KPP: local penetrating shortwave under strong SST damping")
tight_layout()
savefig("local_penetration_mechanism_kpp.png", dpi = 150)
println("saved local_penetration_mechanism_kpp.png")

# -----------------------------------------------------------------------------
# B -- regime map: ΔN²(λ_pen, w₀).  Negative = destratified.
# -----------------------------------------------------------------------------
w0s   = collect(0.0:0.5:3.0) .* 1e-5
λpens = collect(10.0:10.0:60.0)
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
bx.set_title("KPP:  Δ⟨N²⟩ 0–100 m  (heated − no heat)")
colorbar(pc, ax = bx, label = "Δ⟨N²⟩ [s⁻²]")
tight_layout()
savefig("local_penetration_regime_kpp.png", dpi = 150)
println("saved local_penetration_regime_kpp.png")

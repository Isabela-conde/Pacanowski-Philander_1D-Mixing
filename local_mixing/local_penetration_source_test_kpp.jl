# local_penetration_source_test_kpp.jl
# -----------------------------------------------------------------------------
# KPP test that separates the LOCAL route (penetrating shortwave) from the
# REMOTE route (warming the upwelled SOURCE water by ΔT_src), and combines them.
#
# For each regime (FAVOURABLE and REALISTIC) four runs are compared against the
# no-heating baseline:
#   SW only        : penetrating shortwave, source unchanged      (local)
#   source only    : NO surface heat, source warmed by ΔT_src     (remote)
#   SW + source    : both                                         (combined)
#
# Source warming is imposed by restoring the water below z_restore toward
# Tprofile + ΔT_src (ramped in, so no step), representing warm thermocline water
# delivered by the circulation. Same KPP machinery as local_penetration_test_kpp.jl.
# -----------------------------------------------------------------------------

using OceanTurb
const KPP = OceanTurb.KPP
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

z_restore = -100.0        # restore the SOURCE region (deep) toward Tprofile + ΔT_src
d_ramp    = 20.0          # ramp depth for the source-warming anomaly [m]
τ_src     = 60day
z_sub     = -50.0
z_source  = -120.0        # "source water": the thermocline water that gets upwelled

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
# one KPP run.  Q_sfc: penetrating SW (e-folding λ_pen).  ΔT_src: source warming.
# -----------------------------------------------------------------------------
function run_case(; Q_sfc = 0.0, λ_pen = 30.0, λ_fb = 50.0, w₀ = 1.0e-5, ΔT_src = 0.0,
                    N = 128, H = 300.0, Δt = 20minute, tfinal = 365day,
                    wind_stress = 0.01, rec_every = 72)

    grid = UniformGrid(N, H)
    zc   = grid.zc
    Δz   = H / N
    SST₀ = Tprofile(zc[N])
    wz   = [w_prof(z, w₀) for z in zc]
    rad  = Q_sfc > 0 ? [(Q_sfc / (constants.ρ₀ * cp)) * exp(z / λ_pen) / λ_pen for z in zc] :
                       zeros(length(zc))
    # source-restoring target: Tprofile + ΔT_src, ramped from 0 at z_restore
    Ttgt = [Tprofile(z) + ΔT_src * clamp((z_restore - z) / d_ramp, 0.0, 1.0) for z in zc]

    function forcing_T(m, i)
        f = rad[i]
        i > 1              && (f += -wz[i] * (m.solution.T[i] - m.solution.T[i-1]) / Δz)
        zc[i] <= z_restore && (f += (Ttgt[i] - m.solution.T[i]) / τ_src)
        return f
    end

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
    t = Float64[]; sst = Float64[]; Tsub = Float64[]; Tsrc = Float64[]; n2i = Float64[]
    rec!() = (push!(t, model.clock.time); push!(sst, model.solution.T[N]);
              push!(Tsub, T_at(model, z_sub)); push!(Tsrc, T_at(model, z_source));
              push!(n2i, meanN2(model, -100.0, 0.0)))

    rec!()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)
        n % rec_every == 0 && rec!()
    end
    rec!()
    Tfinal = [model.solution.T[i] for i in 1:N]
    return (; t, sst, Tsub, Tsrc, n2i, zc, Tfinal, n2ss = n2i[end])
end

# -----------------------------------------------------------------------------
# figure: three heating scenarios (SW only / source only / both) vs no heating
# -----------------------------------------------------------------------------
function scenario_figure(base, sw, src, both; scheme, λpen, w0, λfb, Q, dTsrc, fname)
    fig, ax = subplots(2, 2, figsize = (13, 9))
    runs = ((sw,   "C3", "SW only (local)"),
            (src,  "C2", @sprintf("source +%.1f°C only (remote)", dTsrc)),
            (both, "C1", "SW + source"))
    zmid, nb = n2_profile(base.Tfinal, base.zc)

    for (r, col, lab) in runs
        ax[1,1].plot(r.t ./ day, r.sst  .- base.sst,  col, label = lab)
        ax[1,2].plot(r.t ./ day, r.Tsrc .- base.Tsrc, col, label = lab)
        ax[2,1].plot(r.Tfinal .- base.Tfinal, r.zc, col, label = lab)
        _, nr = n2_profile(r.Tfinal, r.zc)
        ax[2,2].plot(nr .- nb, zmid, col, label = lab)
    end

    ax[1,1].set_xlabel("time [days]"); ax[1,1].set_ylabel("ΔSST [°C]")
    ax[1,1].set_title("Surface (SST)"); ax[1,1].legend(fontsize = 8)
    ax[1,2].set_xlabel("time [days]"); ax[1,2].set_ylabel("ΔT [°C]")
    ax[1,2].set_title(@sprintf("Source water (%d m)", Int(-z_source))); ax[1,2].legend(fontsize = 8)
    ax[2,1].axvline(0, color = "grey", lw = 0.6); ax[2,1].set_ylim(-150, 0)
    ax[2,1].set_xlabel("ΔT [°C]"); ax[2,1].set_ylabel("z [m]"); ax[2,1].set_title("Final ΔT (vs no heat)")
    ax[2,2].axvline(0, color = "grey", lw = 0.6); ax[2,2].set_ylim(-150, 0)
    ax[2,2].set_xlabel("ΔN² [s⁻²]"); ax[2,2].set_title("Final ΔN² (vs no heat)")

    suptitle(@sprintf("%s:  heating scenarios vs no heating\nλ=%.0f m,  w₀=%.2f m/day,  λ_fb=%.0f W m⁻² K⁻¹,  Q=%.0f W m⁻²,  source +%.1f°C,  2 yr",
                      scheme, λpen, w0*86400, λfb, Q, dTsrc))
    tight_layout(); savefig(fname, dpi = 150); println("saved $fname")
end

function run_regime(; scheme, w0, λfb, λpen, Q, dTsrc, fname)
    base = run_case(Q_sfc = 0.0, w₀ = w0, λ_fb = λfb,                              ΔT_src = 0.0,   Δt = 15minute, tfinal = 730day, rec_every = 96)
    sw   = run_case(Q_sfc = Q,   w₀ = w0, λ_fb = λfb, λ_pen = λpen,                ΔT_src = 0.0,   Δt = 15minute, tfinal = 730day, rec_every = 96)
    src  = run_case(Q_sfc = 0.0, w₀ = w0, λ_fb = λfb,                              ΔT_src = dTsrc, Δt = 15minute, tfinal = 730day, rec_every = 96)
    both = run_case(Q_sfc = Q,   w₀ = w0, λ_fb = λfb, λ_pen = λpen,                ΔT_src = dTsrc, Δt = 15minute, tfinal = 730day, rec_every = 96)
    @printf("%s:  ΔSST  SWonly=%+.2f  srconly=%+.2f  both=%+.2f °C   |   ΔT_source(%dm)  SWonly=%+.2f  srconly=%+.2f  both=%+.2f °C\n",
            scheme, sw.sst[end]-base.sst[end], src.sst[end]-base.sst[end], both.sst[end]-base.sst[end],
            Int(-z_source), sw.Tsrc[end]-base.Tsrc[end], src.Tsrc[end]-base.Tsrc[end], both.Tsrc[end]-base.Tsrc[end])
    scenario_figure(base, sw, src, both; scheme = scheme, λpen = λpen, w0 = w0, λfb = λfb, Q = Q, dTsrc = dTsrc, fname = fname)
end

dTsrc = 1.5

# FAVOURABLE: deep SW, weak upwelling, strong damping
run_regime(scheme = "KPP · favourable", w0 = 0.5e-5,  λfb = 50.0, λpen = 50.0, Q = 80.0, dTsrc = dTsrc,
           fname = "local_source_kpp_favourable.png")

# REALISTIC: shallow SW (~20 m), cold-tongue upwelling (~1.5 m/day), moderate damping
run_regime(scheme = "KPP · realistic", w0 = 1.74e-5, λfb = 30.0, λpen = 20.0, Q = 80.0, dTsrc = dTsrc,
           fname = "local_source_kpp_realistic.png")

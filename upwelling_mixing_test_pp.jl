# upwelling_mixing_test_pp.jl
# -----------------------------------------------------------------------------
# Adds PRESCRIBED UPWELLING to the 1D PP column so we can test the last link in
# the eastern-Pacific argument -- "warmer water was upwelled" -- which the pure
# mixing column (warming_mixing_test_pp.jl) cannot represent.
#
# Same tamed setup (equatorial, f = 0, background-stratified thermocline, gentle
# EUC+SEC shear, PP mixing, SST feedback). NEW ingredients, both operator-split:
#   (1) upwelling w(z) > 0  -> upwind vertical advection of T  (∂T/∂t += −w ∂zT)
#   (2) a deep "source-water" restoring below z_restore toward Tprofile + ΔT_src,
#       representing the remotely-supplied thermocline water. ΔT_src is the knob
#       for "the water that gets upwelled is warmer" (set non-locally in reality).
#
# The point this script makes concrete:
#   * SST in the east is set by a balance between surface heating and upwelling
#     of cold subsurface water -> stronger w gives a colder SST (a "cold tongue").
#   * Warming SST has TWO routes: local surface heating (Q_sfc), and REMOTE
#     source-water warming (ΔT_src) carried up by w. The second route's SST
#     impact scales with w and VANISHES when w = 0 -> it is not a local process.
#
# NOTE: only T is advected by w (the tracer of interest); momentum keeps its
# wind-forced, PP-mixed evolution. That is the standard 1D heat-budget treatment
# and is enough for the SST argument.
# -----------------------------------------------------------------------------

using OceanTurb
using OceanTurb.PacanowskiPhilander: Parameters, KU, KT, local_richardson
using Printf
using PyPlot

@use_pyplot_utils

# ---- physical constants (equator, f = 0) -----------------------------------
constants = Constants(f = 0.0, α = 2.0e-4, β = 8.0e-4, ρ₀ = 1025.0, g = 9.81)
cp = 3991.0   # seawater heat capacity [J kg⁻¹ K⁻¹]

# ---- thermocline + background stratification (same as the tamed mixing test)-
Tdeep, h_th, δ_th = 14.0, 50.0, 30.0
N2_bg = 3.0e-5
Γ_bg  = N2_bg / (constants.g * constants.α)
Tprofile(z, ΔT) = Tdeep + 0.5ΔT * (1 + tanh((z + h_th) / δ_th)) + Γ_bg * z
S₀(z) = 35.0

# ---- fixed two-jet shear: westward SEC over an eastward EUC -----------------
h_core       = 110.0
U_sec, h_sec = -0.25, 50.0
U_euc, σ_euc = 0.60, 45.0
U₀(z) = U_sec * exp(z / h_sec) + U_euc * exp(-((-z - h_core)^2) / (2σ_euc^2))

# ---- upwelling profile: 0 at the surface (divergence), -> w₀ below ----------
# w(z) = w₀ (1 − eᶻ/δw);  z = 0 -> 0,  z = −45 m -> ~0.95 w₀,  deep -> w₀.
δw = 15.0
w_prof(z, w₀) = w₀ * (1 - exp(z / δw))

# ---- deep source-water restoring (the remotely-supplied thermocline) --------
z_restore = -150.0        # restore below this depth
τ_src     = 20day         # restoring timescale

# ---- misc -------------------------------------------------------------------
Ri_crit = 0.25
z_src   = -70.0           # "what gets upwelled": upper-thermocline water

# temperature at the cell centre nearest a target depth
function T_at(model, ztarget)
    zc = nodes(model.solution.T)
    i  = argmin(abs.(zc .- ztarget))
    return model.solution.T[i]
end

# convective adjustment: remove static instabilities (shallower cell colder than
# the one below) by mixing pairs to convergence. Keeps N² ≥ 0 -> Ri ≥ 0, so PP's
# KT = κ₀ + κ₁/(1+5Ri)³ stays positive (it goes NEGATIVE for Ri < −0.2 -> blow-up).
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
    return nothing
end

# -----------------------------------------------------------------------------
# one experiment
#   w₀      upwelling rate [m/s] (+ve upward)
#   Q_sfc   extra surface heating [W/m²]  (local warming route)
#   ΔT_src  source-water warming [°C]     (remote warming route, carried up by w)
# -----------------------------------------------------------------------------
function run_case(; label, w₀ = 0.0, Q_sfc = 0.0, ΔT_src = 0.0, λ_fb = 25.0, λ_pen = 0.0,
                    N = 128, H = 300.0, Δt = 10minute, tfinal = 200day,
                    wind_stress = 0.01, save_profiles = false, save_every = 288)

    # surface heating enters either as a surface flux (λ_pen = 0; only reaches depth
    # by mixing) or as penetrating shortwave ∝ (1/λ_pen)·exp(z/λ_pen) that deposits
    # Q_sfc over an e-folding depth λ_pen (λ_pen > 0). The latter tests whether
    # surface heat can PENETRATE locally to the upwelled-water depth.
    Qᵀ = λ_pen > 0 ? 0.0 : -Q_sfc / (constants.ρ₀ * cp)
    Qᵘ =  wind_stress / constants.ρ₀

    bcs = PacanowskiPhilander.BoundaryConditions(
        FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵘ)),  # U
        ZeroFluxBoundaryConditions(),                                                    # V
        FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵀ)),  # T
        ZeroFluxBoundaryConditions(),                                                    # S
    )

    model = PacanowskiPhilander.Model(grid = UniformGrid(N = N, H = H),
                                      constants = constants,
                                      parameters = Parameters(),
                                      stepper = :BackwardEuler,
                                      bcs = bcs)

    model.solution.U = z -> U₀(z)
    model.solution.V = z -> 0.0
    model.solution.T = z -> Tprofile(z, 10.0)   # ΔT = 10 base; source warming enters via restoring
    model.solution.S = S₀

    Δz     = H / N
    zc     = nodes(model.solution.T)
    SST₀   = model.solution.T[N]
    wz     = [w_prof(z, w₀) for z in zc]                         # upwelling at cell centres
    # penetrating shortwave heating increment per step [°C]; ∫ over depth = Q_sfc
    rad    = λ_pen > 0 ? [(Q_sfc / (constants.ρ₀ * cp)) * exp(z / λ_pen) / λ_pen * Δt for z in zc] :
                         zeros(length(zc))
    # source-restoring target: warming anomaly ramps from 0 at z_restore to full
    # ΔT_src by z_restore − d_ramp, so there is no step inversion at the edge.
    d_ramp = 40.0
    srcanom(z) = ΔT_src * clamp((z_restore - z) / d_ramp, 0.0, 1.0)
    Ttgt   = [Tprofile(z, 10.0) + srcanom(z) for z in zc]

    nsteps = Int(round(tfinal / Δt))
    t = Float64[]; sst = Float64[]; Tsrc = Float64[]
    tsnaps = Float64[]; Tsn = Vector{Float64}[]

    record!() = (push!(t, model.clock.time);
                 push!(sst, model.solution.T[N]);
                 push!(Tsrc, T_at(model, z_src)))
    function snap!()
        save_profiles || return
        push!(tsnaps, model.clock.time)
        push!(Tsn, [model.solution.T[i] for i in 1:N])
    end

    record!(); snap!()
    for n in 1:nsteps
        run_until!(model, Δt, n * Δt)

        # (1) upwind vertical advection of T by w (w>0 upward; upstream = below)
        @inbounds for i in N:-1:2
            model.solution.T[i] += -wz[i] * (model.solution.T[i] - model.solution.T[i-1]) / Δz * Δt
        end
        # (1b) penetrating shortwave heating (λ_pen > 0)
        if λ_pen > 0
            @inbounds for i in 1:N
                model.solution.T[i] += rad[i]
            end
        end
        # (2) deep source-water restoring (maintains / warms the upwelled water)
        @inbounds for i in 1:N
            zc[i] <= z_restore &&
                (model.solution.T[i] += (Δt / τ_src) * (Ttgt[i] - model.solution.T[i]))
        end
        # (3) surface heat-flux feedback (caps SST; net loss as SST rises)
        model.solution.T[N] -= λ_fb * (model.solution.T[N] - SST₀) * Δt /
                               (constants.ρ₀ * cp * Δz)
        # (4) convective adjustment -> no N² < 0, so PP's KT never goes negative
        convect!(model.solution.T, N)

        record!()
        n % save_every == 0 && snap!()
    end

    sstss = sum(sst[t .>= tfinal - 20day]) / count(t .>= tfinal - 20day)   # steady SST
    @printf("[%-8s] w₀=%.1f m/day  Q_sfc=%+5.1f W/m²  ΔT_src=%+.1f°C  ->  SST %.2f->%.2f (steady %.2f) °C\n",
            label, w₀ * 86400, Q_sfc, ΔT_src, sst[1], sst[end], sstss)

    return (; label, t, sst, Tsrc, zc, wz,
              Tfinal = [model.solution.T[i] for i in 1:N], sstss,
              tsnaps, Tsn)
end

# -----------------------------------------------------------------------------
# MAIN RUNS: isolate the two warming routes against a cold-tongue background
# -----------------------------------------------------------------------------
w0   = 2.0e-5     # ~1.7 m/day upwelling
ΔQ   = 80.0       # local surface warming [W/m²]
ΔTsc = 3.0        # remote source-water warming [°C]

noup  = run_case(label = "NOUP",     w₀ = 0.0, save_profiles = true)
upw   = run_case(label = "UPW",      w₀ = w0,  save_profiles = true)
upsfc = run_case(label = "UPW+SFC",  w₀ = w0,  Q_sfc = ΔQ,                save_profiles = true)
upsrc = run_case(label = "UPW+SRC",  w₀ = w0,  ΔT_src = ΔTsc,             save_profiles = true)
upboth= run_case(label = "UPW+BOTH", w₀ = w0,  Q_sfc = ΔQ, ΔT_src = ΔTsc, save_profiles = true)
mains = (noup, upw, upsfc, upsrc, upboth)
cols  = Dict("NOUP"=>"k", "UPW"=>"C0", "UPW+SFC"=>"C1", "UPW+SRC"=>"C3", "UPW+BOTH"=>"C2")

# -----------------------------------------------------------------------------
# FIGURE 1 -- mechanism:  SST(t), final profiles, w(z), and the warm anomaly
#            rising from depth in UPW+SRC
# -----------------------------------------------------------------------------
fig, ax = subplots(2, 2, figsize = (13, 9))

for c in mains
    ax[1,1].plot(c.t ./ day, c.sst, cols[c.label], label = c.label)
end
ax[1,1].axhline(noup.sstss, color = "grey", ls = ":", lw = 1)
ax[1,1].set_xlabel("time [days]"); ax[1,1].set_ylabel("SST [°C]")
ax[1,1].set_title("SST:  upwelling cools, both warming routes warm"); ax[1,1].legend(fontsize = 8)

ax[1,2].plot([Tprofile(z, 10.0) for z in noup.zc], noup.zc, "grey", ls = "--", lw = 1, label = "initial")
for c in mains
    ax[1,2].plot(c.Tfinal, c.zc, cols[c.label], label = c.label)
end
ax[1,2].set_ylim(-250, 0); ax[1,2].set_xlabel("T [°C]"); ax[1,2].set_ylabel("z [m]")
ax[1,2].set_title("final T(z): upwelling lifts isotherms; +SRC warms them"); ax[1,2].legend(fontsize = 7)

# warm anomaly (T − initial) for UPW+SRC and UPW+BOTH, shared colour scale.
# SRC: source heat carried up from depth by w. BOTH: + surface heating up top.
Tinit = [Tprofile(z, 10.0) for z in upw.zc]
Asrc  = reduce(hcat, [snap .- Tinit for snap in upsrc.Tsn])
Aboth = reduce(hcat, [snap .- Tinit for snap in upboth.Tsn])
amax  = max(maximum(abs.(Asrc)), maximum(abs.(Aboth)))
for (a, c, A, ttl) in ((ax[2,1], upsrc,  Asrc,  "UPW+SRC: anomaly T−T₀"),
                       (ax[2,2], upboth, Aboth, "UPW+BOTH: anomaly T−T₀"))
    pc = a.pcolormesh(c.tsnaps ./ day, c.zc, A, cmap = "RdBu_r",
                      vmin = -amax, vmax = amax, shading = "auto")
    a.set_ylim(-250, 0); a.set_xlabel("time [days]"); a.set_ylabel("z [m]")
    a.set_title(ttl); colorbar(pc, ax = a, label = "ΔT [°C]")
end

suptitle("Upwelling closes the argument: the warmer upwelled water must be delivered to the surface by w")
tight_layout()
savefig("upwelling_mechanism.png", dpi = 150)
println("saved upwelling_mechanism.png")

# -----------------------------------------------------------------------------
# FIGURE 2 -- the "not local" plot.  Sweep w four ways: no warming, source
# warming (remote), surface warming (local), and both.  The remote curve peels
# away from the control ONLY as w grows (ΔSST -> 0 at w = 0); the local curve is
# offset even at w = 0.  (Δt = 10 min to match the main runs; NaN runs dropped.)
# -----------------------------------------------------------------------------
w0s = collect(0.0:0.5:4.0) .* 1e-5      # m/s

sst_base = Float64[]; sst_src = Float64[]; sst_sfc = Float64[]; sst_both = Float64[]
for w in w0s
    b  = run_case(label = "sw-base", w₀ = w,                            Δt = 10minute, tfinal = 200day)
    s  = run_case(label = "sw-src",  w₀ = w, ΔT_src = ΔTsc,             Δt = 10minute, tfinal = 200day)
    f  = run_case(label = "sw-sfc",  w₀ = w, Q_sfc  = ΔQ,               Δt = 10minute, tfinal = 200day)
    bo = run_case(label = "sw-both", w₀ = w, Q_sfc = ΔQ, ΔT_src = ΔTsc, Δt = 10minute, tfinal = 200day)
    push!(sst_base, b.sstss); push!(sst_src, s.sstss)
    push!(sst_sfc, f.sstss);  push!(sst_both, bo.sstss)
end
wday  = w0s .* 86400
dsrc  = sst_src  .- sst_base     # ΔSST from remote source warming
dsfc  = sst_sfc  .- sst_base     # ΔSST from local surface warming
dboth = sst_both .- sst_base     # ΔSST from both together
ok    = isfinite.(sst_base) .& isfinite.(sst_src) .& isfinite.(sst_sfc) .& isfinite.(sst_both)

println("\n w₀[m/day]  SST_base  +src   +sfc  +both    Δsrc   Δsfc  Δboth")
for k in eachindex(w0s)
    @printf("  %6.2f   %6.2f  %6.2f %6.2f %6.2f   %+5.2f %+5.2f %+5.2f\n",
            wday[k], sst_base[k], sst_src[k], sst_sfc[k], sst_both[k], dsrc[k], dsfc[k], dboth[k])
end

fig2, bx = subplots(1, 2, figsize = (13, 5))

bx[1].plot(wday[ok], sst_base[ok], "ko-",  label = "no warming")
bx[1].plot(wday[ok], sst_src[ok],  "C3o-", label = @sprintf("source +%.0f°C", ΔTsc))
bx[1].plot(wday[ok], sst_sfc[ok],  "C1o-", label = @sprintf("surface +%.0f W/m²", ΔQ))
bx[1].plot(wday[ok], sst_both[ok], "C2o-", label = "both")
bx[1].set_xlabel("upwelling w₀ [m/day]"); bx[1].set_ylabel("steady SST [°C]")
bx[1].set_title("Steady SST vs upwelling, by warming route"); bx[1].legend(fontsize = 8)

bx[2].plot(wday[ok], dsrc[ok],  "C3o-", label = "source (remote)")
bx[2].plot(wday[ok], dsfc[ok],  "C1o-", label = "surface (local)")
bx[2].plot(wday[ok], dboth[ok], "C2o-", label = "both")
bx[2].axhline(0, color = "grey", lw = 0.7)
bx[2].set_xlabel("upwelling w₀ [m/day]"); bx[2].set_ylabel("ΔSST from warming [°C]")
bx[2].set_title("Remote warming needs w (→ 0 at w = 0); local does not"); bx[2].legend(fontsize = 8)

tight_layout()
savefig("upwelling_sst_sensitivity.png", dpi = 150)
println("saved upwelling_sst_sensitivity.png")

# ------------------------------------------------------------------------------
# FIGURE 3 -- warming-induced anomaly Hovmöllers (T − T_UPW), same colour style
# as the mechanism panel.  Referenced to UPW so the common upwelling cooling is
# removed and the WARMING signal stands out: surface heat stays trapped near the
# surface; source heat is carried UP from depth by w; both combines them.
# ------------------------------------------------------------------------------
anomU(c) = reduce(hcat, [c.Tsn[i] .- upw.Tsn[i] for i in eachindex(c.Tsn)])
hovs  = ((upsfc,  "UPW+SFC − UPW  (local)"),
         (upsrc,  "UPW+SRC − UPW  (remote)"),
         (upboth, "UPW+BOTH − UPW"))
amax3 = maximum(maximum(abs.(anomU(c))) for (c, _) in hovs)

fig3, cx = subplots(1, 3, figsize = (16, 5), sharey = true)
for (k, (c, ttl)) in enumerate(hovs)
    pc = cx[k].pcolormesh(c.tsnaps ./ day, c.zc, anomU(c), cmap = "RdBu_r",
                          vmin = -amax3, vmax = amax3, shading = "auto")
    cx[k].set_ylim(-250, 0); cx[k].set_xlabel("time [days]"); cx[k].set_title(ttl)
    k == 1 && cx[k].set_ylabel("z [m]")
    colorbar(pc, ax = cx[k], label = "ΔT [°C]")
end
suptitle("Warming-induced anomaly (vs UPW): surface heat trapped up top vs source heat carried up by w")
tight_layout()
savefig("upwelling_anomaly_hovmoller.png", dpi = 150)
println("saved upwelling_anomaly_hovmoller.png")

# ------------------------------------------------------------------------------
# FIGURE 4 -- compare experiments: T and N² profiles, all runs overlaid, at the
# start / middle / end (same style as the warming-scenario profiles figure).
# SST is the surface value of T(z); N² = g·α·∂zT (S is constant here).
# ------------------------------------------------------------------------------
function n2_profile(Tcol, zc)
    n = length(zc)
    zmid = zc[2:n-1]
    N2 = [constants.g * constants.α * (Tcol[i+1] - Tcol[i-1]) / (zc[i+1] - zc[i-1])
          for i in 2:n-1]
    return zmid, N2
end

fig4, ax4 = subplots(3, 2, figsize = (11, 14), sharey = true)

im = cld(length(mains[1].tsnaps), 2)     # middle snapshot index
for c in mains
    col = cols[c.label]
    ax4[1,1].plot(c.Tsn[1],   c.zc, col, label = c.label)
    z0, n0 = n2_profile(c.Tsn[1],   c.zc); ax4[1,2].plot(n0, z0, col, label = c.label)
    ax4[2,1].plot(c.Tsn[im],  c.zc, col, label = c.label)
    zM, nM = n2_profile(c.Tsn[im],  c.zc); ax4[2,2].plot(nM, zM, col, label = c.label)
    ax4[3,1].plot(c.Tsn[end], c.zc, col, label = c.label)
    zE, nE = n2_profile(c.Tsn[end], c.zc); ax4[3,2].plot(nE, zE, col, label = c.label)
end

tM = mains[1].tsnaps[im]  / day
tE = mains[1].tsnaps[end] / day
ax4[1,1].set_title("initial:  T(z)   (top = SST)");             ax4[1,2].set_title("initial:  N²(z)")
ax4[2,1].set_title(@sprintf("middle (t = %.0f d):  T(z)", tM)); ax4[2,2].set_title(@sprintf("middle (t = %.0f d):  N²(z)", tM))
ax4[3,1].set_title(@sprintf("final (t = %.0f d):  T(z)",  tE)); ax4[3,2].set_title(@sprintf("final (t = %.0f d):  N²(z)",  tE))
for a in (ax4[1,1], ax4[2,1], ax4[3,1]); a.set_ylabel("z [m]"); a.set_ylim(-250, 0); end
for a in (ax4[1,2], ax4[2,2], ax4[3,2]); a.axvline(0, color = "grey", lw = 0.6); end
ax4[3,1].set_xlabel("T [°C]"); ax4[3,2].set_xlabel("N² [s⁻²]")
for a in (ax4[1,1], ax4[1,2], ax4[2,1], ax4[2,2], ax4[3,1], ax4[3,2]); a.legend(fontsize = 7); end
suptitle("Experiments compared: initial / middle / final  T and N² profiles")
tight_layout()
savefig("upwelling_profiles.png", dpi = 150)
println("saved upwelling_profiles.png")

# ------------------------------------------------------------------------------
# FIGURE 5 -- STABILITY test.  Tests "increased surface warming makes it LESS
# stable" against "warmer source waters change the stability", on stability terms:
#   (a) column stability index  ⟨N²⟩₀₋₁₀₀ₘ  vs time  -> does SFC go up or down?
#   (b,c) ΔN²(z,t) = N²(exp) − N²(UPW) for SFC and SRC -> where/when stability shifts
# This isolates the LOCAL response (1D: surface heat cannot reach the thermocline
# except by mixing), so it shows how much of the claim needs the circulation.
# ------------------------------------------------------------------------------
function meanN2(Tcol, zc, zlo, zhi)
    zmid, N2 = n2_profile(Tcol, zc)
    vals = [N2[i] for i in eachindex(N2) if zlo <= zmid[i] <= zhi]
    return isempty(vals) ? NaN : sum(vals) / length(vals)
end
n2_hov(c) = (n2_profile(c.Tsn[1], c.zc)[1],
             reduce(hcat, [n2_profile(snap, c.zc)[2] for snap in c.Tsn]))

fig5, ex = subplots(1, 3, figsize = (17, 5))

# (a) column stability index vs time
for c in mains
    idx = [meanN2(snap, c.zc, -100.0, 0.0) for snap in c.Tsn]
    ex[1].plot(c.tsnaps ./ day, idx, cols[c.label], label = c.label)
end
ex[1].set_xlabel("time [days]"); ex[1].set_ylabel("⟨N²⟩ 0–100 m [s⁻²]")
ex[1].set_title("Column stability vs time"); ex[1].legend(fontsize = 7)

# (b,c) ΔN²(z,t) vs UPW for SFC and SRC, shared colour scale
zmid, Mupw = n2_hov(upw)
Dsfc = n2_hov(upsfc)[2] .- Mupw
Dsrc = n2_hov(upsrc)[2] .- Mupw
dmax = max(maximum(abs.(Dsfc)), maximum(abs.(Dsrc)))
for (a, D, ttl) in ((ex[2], Dsfc, "ΔN²  SFC − UPW  (surface warming)"),
                    (ex[3], Dsrc, "ΔN²  SRC − UPW  (source warming)"))
    pc = a.pcolormesh(upw.tsnaps ./ day, zmid, D, cmap = "RdBu_r",
                      vmin = -dmax, vmax = dmax, shading = "auto")
    a.set_ylim(-150, 0); a.set_xlabel("time [days]"); a.set_ylabel("z [m]")
    a.set_title(ttl); colorbar(pc, ax = a, label = "ΔN² [s⁻²]")
end
suptitle("Stability response: surface warming vs source warming (relative to UPW).  Red = more stable, blue = less stable")
tight_layout()
savefig("upwelling_stability.png", dpi = 150)
println("saved upwelling_stability.png")

# ------------------------------------------------------------------------------
# FIGURE 6 -- PENETRATION test of the speaker's claim.  Surface heating applied
# as penetrating shortwave with e-folding λ_pen (10/25/50 m): does it reach the
# upwelling source depth and warm the upwelled water, LOCALLY?  Benchmarked
# against SRC, where the source water is warmed directly.
#   (a) final warming anomaly ΔT(z) = T − T_UPW  -> how deep surface heat gets
#   (b) warming at the source depth z_src vs time -> does the upwelled water warm?
# ------------------------------------------------------------------------------
λpens  = (10.0, 25.0, 50.0)
sfcpen = [run_case(label = "SFC λ=$(Int(λ))m", w₀ = w0, Q_sfc = ΔQ, λ_pen = λ,
                   save_profiles = true) for λ in λpens]
pcols  = ("C0", "C1", "C4")

fig6, gx = subplots(1, 2, figsize = (13, 5.5))

# (a) final warming-anomaly profiles (vs UPW)
for (c, pc) in zip(sfcpen, pcols)
    gx[1].plot(c.Tfinal .- upw.Tfinal, c.zc, pc, label = c.label)
end
gx[1].plot(upsrc.Tfinal .- upw.Tfinal, upsrc.zc, "C3", label = @sprintf("SRC (+%.0f°C source)", ΔTsc))
gx[1].axhline(z_src, color = "grey", ls = ":", lw = 1)
gx[1].text(0.03, z_src + 3, "upwelling source", transform = gx[1].get_yaxis_transform(),
           fontsize = 7, color = "grey")
gx[1].set_ylim(-150, 0); gx[1].set_xlabel("ΔT = T − T_UPW  [°C]"); gx[1].set_ylabel("z [m]")
gx[1].set_title("Final warming anomaly: how deep does surface heat get?"); gx[1].legend(fontsize = 7)

# (b) warming reaching the upwelling source depth vs time
for (c, pc) in zip(sfcpen, pcols)
    gx[2].plot(c.t ./ day, c.Tsrc .- upw.Tsrc, pc, label = c.label)
end
gx[2].plot(upsrc.t ./ day, upsrc.Tsrc .- upw.Tsrc, "C3", label = "SRC")
gx[2].axhline(0, color = "grey", lw = 0.7)
gx[2].set_xlabel("time [days]"); gx[2].set_ylabel(@sprintf("ΔT at %d m  [°C]", Int(-z_src)))
gx[2].set_title("Warming reaching the upwelled water"); gx[2].legend(fontsize = 7)

suptitle("Can LOCAL surface-heat penetration warm the upwelled water?  (vs SRC = source warmed directly)")
tight_layout()
savefig("upwelling_penetration.png", dpi = 150)
println("saved upwelling_penetration.png")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/isabelaconde/Documents/GitHub/Pacanowski-Philander_1D-Mixing/local_mixing/"
z = np.linspace(0.0, -300.0, 700)   # depth [m], z<=0

# ---- (1) upwelling  w(z) = w0 (1 - e^{z/δw}) ------------------------------
dw = 15.0
def w_prof(z, w0): return w0 * (1 - np.exp(z / dw))
w_fav = w_prof(z, 0.5e-5)  * 86400.0   # m/day  (favourable)
w_rea = w_prof(z, 1.74e-5) * 86400.0   # m/day  (realistic)

# ---- (2) penetrating shortwave shape  (1/λ) e^{z/λ}  [per m] --------------
def sw_shape(z, lam): return (1.0 / lam) * np.exp(z / lam)

# ---- background temperature profile (what the deep is restored toward) ----
Tdeep, h_th, dth = 14.0, 50.0, 30.0
g, alpha = 9.81, 2.0e-4
Gbg = 3.0e-5 / (g * alpha)
def Tprofile(z): return Tdeep + 0.5 * 10.0 * (1 + np.tanh((z + h_th) / dth)) + Gbg * z
Tbg = Tprofile(z)

def Ttarget(z, z_r, d_ramp, dTsrc):
    T = Tprofile(z) + dTsrc * np.clip((z_r - z) / d_ramp, 0.0, 1.0)
    T = T.copy(); T[z > z_r] = np.nan      # restoring acts only where z <= z_r
    return T

# ===========================================================================
# FIGURE A -- the two prescribed forcing shapes (upwelling, shortwave)
# ===========================================================================
fig, ax = plt.subplots(1, 2, figsize=(9.5, 5.2), sharey=True)

ax[0].plot(w_fav, z, "C0", label="favourable  (w₀=0.43 m/day)")
ax[0].plot(w_rea, z, "C1", label="realistic  (w₀=1.5 m/day)")
ax[0].set_xlabel("upwelling  w(z)  [m/day]"); ax[0].set_ylabel("z [m]")
ax[0].set_title("(1) prescribed upwelling\nw(z) = w₀(1 − e$^{z/δ_w}$),  δ_w=15 m")
ax[0].legend(fontsize=7); ax[0].set_ylim(-250, 0); ax[0].grid(alpha=0.3)

ax[1].plot(sw_shape(z, 20.0), z, "C3", label="λ = 20 m  (realistic)")
ax[1].plot(sw_shape(z, 50.0), z, "C2", label="λ = 50 m  (favourable)")
ax[1].set_xlabel("SW heating shape  (1/λ)e$^{z/λ}$  [m$^{-1}$]")
ax[1].set_title("(2) penetrating shortwave\nfraction of Q deposited per metre")
ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT + "kpp_model_profiles.png", dpi=150)
print("saved kpp_model_profiles.png")

# ===========================================================================
# FIGURES B & C -- the restoring deep ocean, one WITHOUT and one WITH heating.
# Each shows the restoring target T* vs the background profile, in the deep.
# ===========================================================================
def restoring_plot(z_r, dTsrc, color, title, fname):
    fig, ax = plt.subplots(figsize=(5.4, 6))
    Ttgt = Ttarget(z, z_r, 20.0, dTsrc)
    ax.plot(Tbg,  z, "0.5", ls="--", lw=1.6, label="background T (free column)")
    ax.plot(Ttgt, z, color, lw=2.6, label="restoring target  T*")
    if dTsrc > 0:
        ax.fill_betweenx(z, Tbg, Ttgt, where=~np.isnan(Ttgt),
                         color=color, alpha=0.25, label=f"+{dTsrc:.1f} °C anomaly")
    ax.axhspan(-300, z_r, color=color, alpha=0.06)
    ax.axhline(z_r, color=color, ls=":", lw=1.2)
    ax.text(14.8, z_r + 4, f"z$_r$ = {int(z_r)} m", color=color, fontsize=8, ha="right")
    ax.text(14.8, -296, "restoring\nactive below z$_r$", color=color, fontsize=8, ha="right", va="bottom")
    ax.axhline(-120, color="grey", lw=0.6)
    ax.text(9.7, -117, "source depth (−120 m)", color="grey", fontsize=7, va="bottom")
    ax.set_xlim(9.5, 15.0); ax.set_ylim(-300, -40)
    ax.set_xlabel("temperature [°C]"); ax.set_ylabel("z [m]")
    ax.set_title(title); ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(OUT + fname, dpi=150); print("saved", fname)

restoring_plot(-200.0, 0.0, "C0",
               "Restoring: cold anchor (no source heating)\nlocal_penetration_test_kpp.jl",
               "kpp_restoring_no_source.png")
restoring_plot(-100.0, 1.5, "C3",
               "Restoring: heated source (+1.5 °C)\nlocal_penetration_source_test_kpp.jl",
               "kpp_restoring_heated_source.png")

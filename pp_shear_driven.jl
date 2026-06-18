# using OceanTurb

# using OceanTurb.PacanowskiPhilander: Parameters

# @use_pyplot_utils # add utilities for plotting OceanTurb Fields


using OceanTurb
 
using OceanTurb.PacanowskiPhilander: Parameters, KU, KT, local_richardson

@use_pyplot_utils
# using PyPlot   # provides subplots, plot, sca, cla, xlabel, legend, L"...", etc.
 
# # Plotting helpers that used to come from @use_pyplot_utils (removed from OceanTurb).
# # Let `plot` accept an OceanTurb Field directly, and add the removespines helper.
# import PyPlot: plot
# plot(f::OceanTurb.AbstractField, args...; kwargs...) =
#     plot(OceanTurb.data(f), OceanTurb.nodes(f), args...; kwargs...)
# removespine(side, ax=gca()) = ax.spines[side].set_visible(false)
# removespines(sides...; ax=gca()) = for side in sides; removespine(side, ax); end
 
function makeplot!(fig, axs, model)
 
    fig.suptitle("\$ t = \$ $(prettytime(model.clock.time))")
 
    # PP diagnostics live at cell faces (they depend on vertical gradients)
    Ri = FaceField(model.grid)
    ν  = FaceField(model.grid)   # eddy viscosity   = KU
    κ  = FaceField(model.grid)   # eddy diffusivity = KT
 
    for i in eachindex(Ri)
        Ri[i] = local_richardson(model, i)
        ν[i]  = KU(model, i)
        κ[i]  = KT(model, i)
    end
 
    markerkwargs = Dict(:marker=>"s", :markersize=>1)
 
    sca(axs[1]); cla()
    plot(model.solution.T; linestyle="-", markerkwargs...)
    removespines("top", "right")
    xlabel("Temperature (\$ {}^\\circ \\mathrm{C} \$)")
    ylabel(L"z \, \mathrm{(m)}")
 
    sca(axs[2]); cla()
    plot(model.solution.U; label=L"U", linestyle="-",  markerkwargs...)
    plot(model.solution.V; label=L"V", linestyle="--", markerkwargs...)
    removespines("top", "right", "left")
    legend()
    xlabel("Velocity (m s\$^{-1}\$)")
 
    sca(axs[3]); cla()
    plot(Ri; linestyle="-", markerkwargs...)
    removespines("top", "right", "left")
    xlabel(L"Ri")
    xlim(-1, 10)   # Ri can blow up where shear → 0; clip for readability
 
    sca(axs[4]); cla()
    plot(ν; label=L"\nu", linestyle="-",  markerkwargs...)
    plot(κ; label=L"\kappa", linestyle="--", markerkwargs...)
    removespines("top", "right", "left")
    legend()
    xlabel("Diffusivity (m\$^2\$ s\$^{-1}\$)")
 
    for i = 2:length(axs)
        axs[i].tick_params(left=false, labelleft=false)
    end
 
    ylim(-H, 1)
    pause(0.1)
 
    return nothing
end
 
constants = Constants(f=1e-4)
 
 N = 128        # Model resolution
 H = 128        # Vertical extent of the model domain
Qᵘ = -1e-4      # Surface momentum flux (wind stress / ρ), drives shear
N² = 1e-5       # Interior/initial buoyancy gradient
Δt = 1minute
 
dTdz = N² / (constants.α * constants.g)
 
parameters = Parameters(
    Cν₀ = 1e-4,   # background viscosity ν₀
    Cν₁ = 1e-2,   # max additional viscosity ν₁
    Cκ₀ = 1e-5,   # background diffusivity κ₀
    Cκ₁ = 1e-2,   # max additional diffusivity κ₁
    Cc  = 5.0,    # Richardson-number coefficient
    Cn  = 2.0,    # exponent n
)
 
bcs = PacanowskiPhilander.BoundaryConditions(
    FieldBoundaryConditions(FluxBoundaryCondition(0.0), FluxBoundaryCondition(Qᵘ)),       # U
    ZeroFluxBoundaryConditions(),                                                          # V
    FieldBoundaryConditions(GradientBoundaryCondition(dTdz), FluxBoundaryCondition(0.0)),  # T
    ZeroFluxBoundaryConditions(),                                                          # S
)

model = PacanowskiPhilander.Model(grid = UniformGrid(N=N, H=H),
                             constants = constants,
                            parameters = parameters,
                               stepper = :ForwardEuler,
                                   bcs = bcs)

T₀(z) = 20 + dTdz * z
model.solution.T = T₀

# Run the model
fig, axs = subplots(ncols=4, figsize=(12, 5), sharey=true)
 
for iplot = 1:12
    run_until!(model, Δt, iplot * 1hour)
    OceanTurb.update!(model)
    makeplot!(fig, axs, model)
end
 
fig.savefig("model_profiles.png", dpi=150, bbox_inches="tight")
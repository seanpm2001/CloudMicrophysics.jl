import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const CMO = CM.Common
const CMI = CM.HetIceNucleation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
# Boiler plate code to have access to model parameters and constants
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
prs = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(prs)

# Knopf and Alpert 2013 Figure 5A: Cirrus
# https://doi.org/10.1039/C3FD00035D 
temp = [228.20357, 228.33571, 228.50357, 228.75000, 228.92143, 229.16786, 229.39643, 229.52143]
KA13_fig5A_J = [70170382.86704, 12101528.74384, 1277935.32665, 52710.05764, 6040.39151, 293.39697, 18.97171, 4.18121]

# Our parameterization
dust_type = CMT.IlliteType()
x_sulph = collect(0.008:0.025:0.2)
x_sulph_const = 0.1

Δa_w = Vector{Float64}(undef, length(temp))
J_ABIFM = Vector{Float64}(undef, length(temp))
constx_Δa_w = Vector{Float64}(undef, length(temp))
constx_J_ABIFM = Vector{Float64}(undef, length(temp))

it = 1
for T in temp
    #Δa_w[it] = CMO.Delta_a_w(prs, x_sulph, T)

    #a_ice = exp((210368 + 131.438*T - (3.32373e6 /T) - 41729.1*log(T))/(8.31441*T)) # exp(large negative) = 0.0 right now
    a_ice = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice()) / CMO.H2SO4_soln_saturation_vapor_pressure(FT(0),T)
    a_sol = CMO.H2SO4_soln_saturation_vapor_pressure(x_sulph[it],T) / CMO.H2SO4_soln_saturation_vapor_pressure(FT(0),T)

    Δa_w[it] = max(a_sol - a_ice, FT(0.))
    println(Δa_w[it])
    J_ABIFM[it] = CMI.ABIFM_J(dust_type, Δa_w[it]) / 100^2

    # plotting w constant x_sulph for comparison
    constx_a_sol = CMO.H2SO4_soln_saturation_vapor_pressure(x_sulph_const,T) / CMO.H2SO4_soln_saturation_vapor_pressure(FT(0),T)
    constx_Δa_w[it] = max(constx_a_sol - a_ice, FT(0.))
    constx_J_ABIFM[it] = CMI.ABIFM_J(dust_type, constx_Δa_w[it]) / 100^2

    global it += 1
end

# Plot results
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "J_het [cm^-2 s^-1]", xlabel = "Temp [K]", yscale = log10)
MK.ylims!(ax1, FT(1), FT(1e8))

MK.lines!(ax1, temp, KA13_fig5A_J, label = "KA13 Figure 5A")
MK.lines!(ax1, temp, J_ABIFM, label = "hackie x_sulph")
MK.lines!(ax1, temp, constx_J_ABIFM, label = "constant x_sulph")
MK.save("ABIFM_cirrus.svg", fig)

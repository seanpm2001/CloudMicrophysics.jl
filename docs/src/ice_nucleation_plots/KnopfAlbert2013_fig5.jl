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
x_sulph = FT(0.1)

Δa_w = Vector{Float64}(undef, length(temp))
J_ABIFM = Vector{Float64}(undef, length(temp))

it = 1
for T in temp
    Δa_w[it] = CMO.Delta_a_w(prs, x_sulph, T)
    J_ABIFM[it] = CMI.ABIFM_J(dust_type, Δa_w[it]) / 100^2
    global it += 1
end

# Plot results
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "J_het [cm^-2 s^-1]", xlabel = "Temp [K]", yscale = log10)
MK.ylims!(ax1, FT(1), FT(1e8))

MK.lines!(ax1, temp, KA13_fig5A_J)
MK.lines!(ax1, temp, J_ABIFM)
MK.save("ABIFM_cirrus.svg", fig)

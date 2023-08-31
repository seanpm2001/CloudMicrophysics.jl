import CairoMakie as MK
import NLsolve as NL
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

# Melting point curve (eq. 3 from Swanson 2008)
function awm_swanson(T)
    return exp(15.8083 + 25301.4 / T - 399752.0 / T / T - 5018.85 * log(T)/T)
end

# water activity of the solution - the ratio between
# the water vapour pressures of the solution
# and of pure water under the same conditions
function water_act_default(x, T)
    p_sol = CMO.H2SO4_soln_saturation_vapor_pressure(prs, x, T)
    p_sat_liq = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    return p_sol / p_sat_liq
end
function water_act_new(x, T)
    p_sol = CMO.H2SO4_soln_saturation_vapor_pressure(prs, x, T)
    p_pure = CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.0, T)
    return p_sol / p_pure
end

# wide temperature range to plot things
T_range = range(180, stop=280, length=100)

# our saturtion pressures over liquid and ice
p_sat_liq = [TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) for T in T_range]
p_sat_ice = [TD.saturation_vapor_pressure(thermo_params, T, TD.Ice()) for T in T_range]
# equilibrium pressure over solution (Luo 1995)
p_sol_50 = [CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.5, T) for T in T_range]
p_sol_30 = [CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.3, T) for T in T_range]
p_sol_10 = [CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.1, T) for T in T_range]
p_sol_0 = [CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.0, T) for T in T_range]

lw = 3
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "p [hPa]", xlabel = "T [C]", limits = ((-90, -35), (0, 0.25)))
#MK.lines!(ax1, T_range, p_sat_cool ./ p_sat_liq, label = "sat ice cool", color = :blue)
MK.lines!(ax1, T_range .- 273.15, p_sat_liq ./ 100, linewidth = lw, label = "p sat liq", color = :blue)
MK.lines!(ax1, T_range .- 273.15, p_sat_ice ./ 100, linewidth = lw, label = "p sat ice", color = :cyan)
MK.lines!(ax1, T_range .- 273.15, p_sol_50 ./ 100,  linewidth = lw, label = "p sol 0.5", color = :purple)
MK.lines!(ax1, T_range .- 273.15, p_sol_30 ./ 100,  linewidth = lw, label = "p sol 0.3", color = :red)
MK.lines!(ax1, T_range .- 273.15, p_sol_10 ./ 100,  linewidth = lw, label = "p sol 0.1", color = :pink)
MK.lines!(ax1, T_range .- 273.15, p_sol_0  ./ 100,  linewidth = lw, label = "p sol 0",   color = :orange)
MK.vlines!(ax1, [185-273.15, 235-273.15], color = :gray)
MK.axislegend(position = :lt )
MK.save("saturation_pressures_liq_ice_corr_range.pdf", fig, pt_per_unit = 2)

# Trying to plot Fig 3 from Swanson
function melt_point_curve(x_sol)
    function melt_point_problem(T)
         x = x_sol
         p_ice_vect = Vector{Float64}(undef, length(T))
         p_liq_vect = Vector{Float64}(undef, length(T))
         p_pur_vect = Vector{Float64}(undef, length(T))
         p_sol_vect = Vector{Float64}(undef, length(T))
         @. p_ice_vect = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
         @. p_liq_vect = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
         @. p_pur_vect = CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.0, T)
         @. p_sol_vect = CMO.H2SO4_soln_saturation_vapor_pressure(prs, x, T)
         #return p_ice_vect .- p_sol_vect
         #return p_liq_vect .- p_sol_vect
         return p_pur_vect .- p_sol_vect
    end
    sol = NL.nlsolve(melt_point_problem, [270.0,])
    return sol.zero[1]
end

#T_range = range(190, stop=240, length=100)
x_range = range(0.1, stop=0.4, length=200)
T_melt = [melt_point_curve(x_sol) for x_sol in x_range]

a_wm_swan = [awm_swanson(T) for T in T_range]
a_wm_swan_trans = a_wm_swan .+ 0.32
a_w_ice = p_sat_ice ./ p_sat_liq
a_w_ice_trans = p_sat_ice ./ p_sat_liq .+ 0.32

#println(x_range)
#println(T_melt)
#println(water_act_new(0.15, 278.1))
#println(water_act_new.(x_range, T_melt))

fig = MK.Figure(resolution = (700, 500))
#ax1 = MK.Axis(fig[1, 1], ylabel = "T [K]", xlabel = "-aw", limits = ((-1, -0.6), (190, 240)))
ax1 = MK.Axis(fig[1, 1], ylabel = "T [K]", xlabel = "-aw", limits = (nothing, nothing))


#MK.lines!(ax1, a_w1 .* -1, T_range, label = "CM default", linestyle = :dash, color = :red)
#MK.lines!(ax1, a_w1_new .* -1, T_range, label = "Using p(0,T)", color = :red)

MK.lines!(ax1, water_act_new.(x_range, T_melt) .* -1, T_melt, label = "trying", color = :red)

MK.lines!(ax1, a_w_ice .* -1,         T_range, label = "a_w_ice",       color = :blue, linestyle =:dash)
MK.lines!(ax1, a_w_ice_trans .* -1,   T_range, label = "a_w_ice_trans", color = :blue)
MK.lines!(ax1, a_wm_swan .* -1,       T_range, label = "a_wm_swanson",  color = :green, linestyle =:dash)
MK.lines!(ax1, a_wm_swan_trans .* -1, T_range, label = "a_wm_swanson_trans",  color = :green)
MK.axislegend()
MK.save("Swanson_Fig3.pdf", fig, pt_per_unit = 2)

#
#
#T_melt = [melt_point_curve(x_sol) for x_sol in x_range]
#fig = MK.Figure(resolution = (800, 600))
#ax1 = MK.Axis(fig[1, 1], ylabel = "T [C]", xlabel = "x [%]")
#MK.lines!(ax1, x_range .* 100, T_melt, color = :blue)
#fig

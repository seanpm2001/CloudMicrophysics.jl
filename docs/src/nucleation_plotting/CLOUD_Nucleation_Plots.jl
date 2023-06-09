using Plots
using CLIMAParameters

include("../../src/Nucleation.jl")
using .Nucleation

# Testing for CLOUD-experiment based nucleation rates.

FT = Float64
toml_dict = CLIMAParameters.create_toml_dict(FT)
param_names = [
    "u_b_n",
    "v_b_n",
    "w_b_n",
    "u_b_i",
    "v_b_i",
    "w_b_i",
    "u_t_n",
    "v_t_n",
    "w_t_n",
    "u_t_i",
    "v_t_i",
    "w_t_i",
    "p_t_n",
    "p_A_n",
    "a_n",
    "p_t_i",
    "p_A_i",
    "a_i",
    "p_b_n",
    "p_b_i",
    "p_t_n",
    "p_t_i",
]
params = CLIMAParameters.get_parameter_values!(toml_dict, param_names)
params = (; params...)

function plot_pure_h2so4_nucleation_rate(
    h2so4_concentrations,
    nh3_conc,
    negative_ion_conc,
    temp,
    params,
)
    rates = map(h2so4_concentrations) do h2so4_conc
        sum(
            Nucleation.h2so4_nucleation_rate(
                h2so4_conc * 1e6,
                nh3_conc,
                negative_ion_conc,
                temp,
                params,
            ),
        ) * 1e-6
    end
    Plots.plot!(
        # title = title,
        h2so4_concentrations,
        rates,
        xaxis = :log,
        yaxis = :log,
        lw = 3,
        ylims = (1e-4, 1e3),
        ylabel = "Nucleation rate (cm⁻³ s⁻¹)",
        xlabel = "[H2SO4] (cm⁻³)",
        label = "$temp K",
        palette = :Dark2_5,
    )
    return rates

end

Plots.plot()
h2so4_concs = [
    10 .^ (5:0.125:7.5),
    10 .^ (5.5:0.125:8),
    10 .^ (5.5:0.125:9),
    10 .^ (7.0:0.125:9.3),
    10 .^ (7:0.125:9.3),
]
temps = [208, 223, 248, 278, 292]
for (temp, h2so4_conc) in zip(temps, h2so4_concs)
    plot_pure_h2so4_nucleation_rate(h2so4_conc, 0, 0, temp, params)
end

dunne_points = (
    [
        (1693569.5886516434, 0.10744256664074078),
        (4555751.638117136, 0.21435379525118536),
        (11230760.811240233, 8.027961753853548),
        (8181937.684310501, 8.9829327833832),
        (4102562.698490346, 7.1081003572394215),
        (4835138.601681659, 1.2603227673494075),
        (2018248.1572122737, 0.04569896206007829),
        (16405303.798002133, 160.21543147865268),
        (7469878.252555515, 26.26346644150873),
        (1496299.43370972, 0.02894783792968937),
        (10450026.858095868, 27.95569539399296),
    ],
    [
        (999424.1083310976, 0.00013298710250629),
        (6687403.049764221, 0.034145488738335936),
        (5209507.122917159, 0.08030857221391505),
        (5661744.37240141, 0.3197141845527293),
        (6641171.525602876, 0.38102404299462694),
        (12748176.303258475, 2.017260128780402),
        (9137630.182549853, 3.2680275894101185),
        (13289996.649555974, 16.200328072641284),
        (28308721.8620071, 95.70891236771261),
    ],
    [
        (8254041.852680173, 0.012270125831294304),
        (13768571.648527596, 0.0177583214837628),
        (18214475.36395946, 0.031408181821309146),
        (14677992.676220706, 0.02363404247446858),
        (34254861.99786931, 0.1819567210692885),
        (32132494.68997057, 1.2939233771410916),
    ],
    [
        (262469829.5872748, 0.022542924011446547),
        (377277947.86104393, 0.038095843909951424),
        (387013239.61175257, 0.09004715173678965),
        (538477507.2128237, 0.2338330660813884),
        (775206678.5767579, 2.2065295605641566),
        (1488686579.5998137, 18.238673966646644),
        (315040901.1713288, 0.1052524521566819),
        (652177614.9057266, 2.58010582090349),
        (504438120.07904917, 0.4503914170005437),
    ],
)

for points in dunne_points
    Plots.plot!(points, label = "", seriestype = :scatter)
end
Plots.plot!(xticks = 10 .^ (5:10))
Plots.svg("CLOUD_nucleation")
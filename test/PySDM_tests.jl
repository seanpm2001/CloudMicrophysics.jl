import Test

import CloudMicrophysics
import CLIMAParameters

const TT = Test
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CP = CLIMAParameters

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

# Atmospheric conditions
T = [282.8054379700115, 282.82011913510763, 282.8182541470769, 282.816547286481,
     282.8152302772072, 282.83225081516434, 282.83107153716537, 282.8300594435075,
     282.82918785057166, 282.8284346224201, 282.8277816799783, 282.8272141498667,
     282.8267193163115, 282.8262868328548, 282.8259080527273, 282.825575586196,
     282.8252831734087, 282.8250255282393, 282.82479813646944, 282.8245971038684]     # air temperature
p = [99352.28143153602, 99381.95622299137, 99381.94656495028, 99381.85770219186,
     99381.90580284428, 99411.64002911752, 99411.58811377818, 99411.54174606412,
     99411.50019030878, 99411.46280238577, 99411.42915075587, 99411.39885950244,
     99411.37147577599, 99411.34671755912, 99411.32435993387, 99411.30416091887,
     99411.28589558789, 99411.26938254002, 99411.25445547458, 99411.24096135839] # air pressure
w = .5        # vertical velocity

r_factors = [0.26666667, 0.29333333, 0.32, 0.34666667, 0.37333333, 0.4,
             0.42666667, 0.45333333, 0.48, 0.50666667, 0.53333333, 0.56,
             0.58666667, 0.61333333, 0.64, 0.66666667, 0.69333333, 0.72,
             0.74666667, 0.77333333]

r_dry = 0.04 * 1e-6
N = 1000. * 1e6
stdev_1 = 1.4
stdev_2 = 1.6

κ_sulfate = .53
M_sulfate = .12314
ρ_sulfate = 1770.0
ϕ_sulfate = 1.0
ν_sulfate = 3.0
ϵ_sulfate = 1.0

mode1_κ = AM.Mode_κ(
    r_dry,
    stdev_1,
    N,
    (1.0,),
    (M_sulfate,),
    (κ_sulfate,),
    1
)
mode1_B = AM.Mode_B(
    r_dry,
    stdev_1,
    N,
    (1.0,),
    (ϵ_sulfate,),
    (ϕ_sulfate,),
    (M_sulfate,),
    (ν_sulfate,),
    (ρ_sulfate,),
    1,
)

N_frac_κ = zeros(length(r_factors))
N_frac_B = zeros(length(r_factors))

for i in 1:length(r_factors)

    mode2_κ = AM.Mode_κ(
        r_dry/r_factors[i],
        stdev_2,
        N,
        (1.0,),
        (M_sulfate,),
        (κ_sulfate,),
        1
    )
    mode2_B = AM.Mode_B(
        r_dry,
        stdev_1,
        N,
        (1.0,),
        (ϵ_sulfate,),
        (ϕ_sulfate,),
        (M_sulfate,),
        (ν_sulfate,),
        (ρ_sulfate,),
        1,
    )

    model_κ = AM.AerosolDistribution((mode1_κ, mode2_κ))
    model_B = AM.AerosolDistribution((mode1_B, mode2_B))

    N_total_κ = model_κ.Modes[1].N + model_κ.Modes[2].N
    N_total_B = model_B.Modes[1].N + model_B.Modes[2].N

    N_act_κ = AA.total_N_activated(param_set, model_κ, T[i], p[i], w)
    N_act_B = AA.total_N_activated(param_set, model_B, T[i], p[i], w)

    N_frac_κ[i] = N_act_κ/N_total_κ
    N_frac_B[i] = N_act_B/N_total_B

end

println("---")
println(N_frac_κ)
println("---")
println(N_frac_B)
println("---")


import Test

import CloudMicrophysics
import CLIMAParameters

const TT = Test
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CP = CLIMAParameters

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

Ts = [282.8054379700115, 282.8193237698699, 282.8172571225181, 282.815027266051,
      282.8316384654537, 282.8301777589934, 282.82898919258645, 282.8280157791203,
      282.8272141498667, 282.82655032735556, 282.8259981385378, 282.82553695943926,
      282.8251503057445, 282.82482504355147, 282.82455056951125, 282.8243182723469,
      282.82412115100965, 282.82395344176706]

ps = [99352.28143153602, 99381.87961549801, 99382.01054047632, 99381.89910729272,
      99411.61334429379, 99411.5472583031, 99411.49045304587, 99411.44135605945,
      99411.39885950244, 99411.36189769715, 99411.32973711056, 99411.30177613927,
      99411.27743121835, 99411.25624008477, 99411.23779596454, 99411.22174453341,
      99411.20778606429, 99411.19564383295]

r_factors = [0.26666667, 0.30333333, 0.34, 0.37666667, 0.41333333, 0.45,
             0.48666667, 0.52333333, 0.56, 0.59666667, 0.63333333, 0.67,
             0.70666667, 0.74333333, 0.78, 0.81666667, 0.85333333, 0.89]
w = .5

κ_as = .53
r_dry = 4 * 1e-8
N_mode = 100000.
s_geom1 = 1.4
s_geom2 = 1.6
M_as = 0.13214

AM1 = AM.Mode_κ(
    r_dry,
    s_geom1,
    N_mode,
    (1.0, ),
    (M_as,),
    (κ_as,),
    1,
)
N_frac = zeros(length(r_factors))
for i in 1:length(r_factors)
    T = Ts[i]
    p = ps[i]
    r_factor = r_factors[i]
    AM2 = AM.Mode_κ(
        r_dry/r_factor,
        s_geom2,
        N_mode,
        (1.0, ),
        (M_as,),
        (κ_as,),
        1,
    )
    Model = AM.AerosolDistribution((AM1, AM2))
    tot_act = AA.total_N_activated(param_set, Model, T, p, w)
    println("Act no.: ", tot_act)
    N_tot = 2*N_mode
    N_frac[i] = tot_act/N_tot
end
println("Fraction Activated")
println(N_frac)

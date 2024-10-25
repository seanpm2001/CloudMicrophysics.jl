"""
    Module for functions shared by different parameterizations.
"""
module Common

import SpecialFunctions as SF

import Thermodynamics as TD
const TPS = TD.Parameters.ThermodynamicsParameters

import ..Parameters as CMP
const HPS = CMP.H2SO4SolutionParameters

export G_func
export H2SO4_soln_saturation_vapor_pressure
export a_w_xT
export a_w_eT
export a_w_ice
export Chen2022_monodisperse_pdf
export Chen2022_exponential_pdf
export Chen2022_vel_coeffs_B1
export Chen2022_vel_coeffs_B2
export Chen2022_vel_coeffs_B4

"""
    G_func(air_props, tps, T, Liquid())
    G_func(air_props, tps, T, Ice())

 - `air_props` - struct with air parameters
 - `tps` - struct with thermodynamics parameters
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Utility function combining thermal conductivity and vapor diffusivity effects.
"""
function G_func(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TPS,
    T::FT,
    ::TD.Liquid,
) where {FT}
    R_v = TD.Parameters.R_v(tps)
    L = TD.latent_heat_vapor(tps, T)
    p_vs = TD.saturation_vapor_pressure(tps, T, TD.Liquid())

    return FT(1) /
           (L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs)
end
function G_func(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TPS,
    T::FT,
    ::TD.Ice,
) where {FT}
    R_v = TD.Parameters.R_v(tps)
    L = TD.latent_heat_sublim(tps, T)
    p_vs = TD.saturation_vapor_pressure(tps, T, TD.Ice())

    return FT(1) /
           (L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs)
end

"""
    A Heaviside step function
"""
heaviside(x::FT) where {FT} = FT(x > 0)

"""
    logistic_function(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the curve, characterizing steepness of the transition

Returns the value of the logistic function for smooth transitioning at thresholds. This is
a normalized curve changing from 0 to 1 while x varies from 0 to Inf (for positive k). For
x < 0 the value at x = 0 (zero) is returned. For x_0 = 0 H(x) is returned.
"""
function logistic_function(x::FT, x_0::FT, k::FT) where {FT}

    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return FT(1)
    end

    return FT(1) / (FT(1) + exp(-k * (x / x_0 - x_0 / x)))
end

"""
    logistic_function_integral(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the logistic function, characterizing steepness of the transition

Returns the value of the indefinite integral of the logistic function, for smooth transitioning
of piecewise linear profiles at thresholds. This curve smoothly transition from y = 0
for 0 < x < x_0 to y = x - x_0 for x_0 < x.
"""
function logistic_function_integral(x::FT, x_0::FT, k::FT) where {FT}
    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return x
    end

    # translation of the curve in x and y to enforce zero at x = 0
    trnslt = -log(FT(1) - exp(-k)) / k

    kt = k * (x / x_0 - FT(1) + trnslt)
    return (kt > FT(40)) ? x - x_0 : (log(FT(1) + exp(kt)) / k - trnslt) * x_0
end

"""
    H2SO4_soln_saturation_vapor_pressure(prs, x, T)

 - `prs` - a struct with H2SO4 solution free parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns the saturation vapor pressure above a sulphuric acid solution droplet in Pa.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight
"""
function H2SO4_soln_saturation_vapor_pressure(
    (;
        T_max,
        T_min,
        w_2,
        c1,
        c2,
        c3,
        c4,
        c5,
        c6,
        c7,
    )::CMP.H2SO4SolutionParameters{FT},
    x::FT,
    T::FT,
) where {FT}

    @assert T < T_max
    @assert T > T_min

    w_h = w_2 * x
    p_sol =
        exp(
            c1 - c2 * x + c3 * x * w_h - c4 * x * w_h^2 +
            (c5 + c6 * x - c7 * x * w_h) / T,
        ) * 100 # * 100 converts mbar --> Pa
    return p_sol
end

"""
    a_w_xT(H2SO4_prs, tps, x, T)

 - `H2SO4_prs` - a struct with H2SO4 solution free parameters
 - `tps` - a struct with thermodynamics parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns water activity of H2SO4 containing droplet.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight.
"""
function a_w_xT(H2SO4_prs::HPS, tps::TPS, x::FT, T::FT) where {FT}

    p_sol = H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x, T)
    p_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())

    return p_sol / p_sat
end

"""
    a_w_eT(tps, e, T)

 - `tps` - struct with thermodynamics parameters
 - `e` - partial pressure of water [Pa]
 - `T` - air temperature [K].

Returns water activity of pure water droplet.
Valid when droplet is in equilibrium with surroundings.
"""
function a_w_eT(tps::TPS, e::FT, T::FT) where {FT}
    # RH
    return e / TD.saturation_vapor_pressure(tps, T, TD.Liquid())
end

"""
    a_w_ice(tps, T)

- `tps` - struct with thermodynamics parameters
 - `T` - air temperature [K].

Returns water activity of ice.
"""
function a_w_ice(tps::TPS, T::FT) where {FT}

    return TD.saturation_vapor_pressure(tps, T, TD.Ice()) /
           TD.saturation_vapor_pressure(tps, T, TD.Liquid())
end

"""
    Chen2022_vel_coeffs_B1(coeffs, ρₐ)

 - coeffs - a struct with terminal velocity free parameters
 - ρ - air density

Returns the coefficients from Table B1 Appendix B in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function Chen2022_vel_coeffs_B1(
    coeffs::CMP.Chen2022VelTypeRain,
    ρₐ::FT,
) where {FT}
    (; ρ0, a, a3_pow, b, b_ρ, c) = coeffs
    # Table B1
    q = exp(ρ0 * ρₐ)
    ai = (a[1] * q, a[2] * q, a[3] * q * ρₐ^a3_pow)
    bi = (b[1] - b_ρ * ρₐ, b[2] - b_ρ * ρₐ, b[3] - b_ρ * ρₐ)
    ci = (c[1], c[2], c[3])
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return (aiu, bi, ciu)
end

"""
    Chen2022_vel_coeffs_B2(coeffs, ρₐ, ρᵢ)

 - coeffs - a struct with terminal velocity free parameters
 - ρₐ - air density
 - ρᵢ - apparent density of ice particles

Returns the coefficients from Table B2 Appendix B in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function Chen2022_vel_coeffs_B2(
    coeffs::CMP.Chen2022VelTypeSmallIce,
    ρₐ::FT,
    ρᵢ::FT,
) where {FT}
    (; A, B, C, E, F, G) = coeffs
    # Table B3
    As = A[2] * (log(ρᵢ))^2 − A[3] * log(ρᵢ) + A[1]
    Bs = FT(1) / (B[1] + B[2] * log(ρᵢ) + B[3] / sqrt(ρᵢ))
    Cs = C[1] + C[2] * exp(C[3] * ρᵢ) + C[4] * sqrt(ρᵢ)
    Es = E[1] - E[2] * (log(ρᵢ))^2 + E[3] * sqrt(ρᵢ)
    Fs = -exp(F[1] - F[2] * (log(ρᵢ))^2 + F[3] * log(ρᵢ))
    Gs = FT(1) / (G[1] + G[2] / (log(ρᵢ)) - G[3] * log(ρᵢ) / ρᵢ)
    # Table B2
    ai = (Es * ρₐ^As, Fs * ρₐ^As)
    bi = (Bs + ρₐ * Cs, Bs + ρₐ * Cs)
    ci = (FT(0), Gs)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return (aiu, bi, ciu)
end

"""
    Chen2022_vel_coeffs_B4(coeffs, ρₐ, ρᵢ)

 - coeffs - a struct with terminal velocity free parameters
 - ρₐ - air density
 - ρᵢ - apparent density of ice particles

Returns the coefficients from Table B4 Appendix B in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function Chen2022_vel_coeffs_B4(
    coeffs::CMP.Chen2022VelTypeLargeIce,
    ρₐ::FT,
    ρᵢ::FT,
) where {FT}
    (; A, B, C, E, F, G, H) = coeffs
    # Table B5
    Al = A[1] + A[2] * log(ρᵢ) + A[3] * (ρᵢ)^(-3 / 2)
    Bl = exp(B[1] + B[2] * log(ρᵢ)^2 + B[3] * log(ρᵢ))
    Cl = exp(C[1] + C[2] / log(ρᵢ) + C[3] / ρᵢ)
    El = E[1] + E[2] * log(ρᵢ) * sqrt(ρᵢ) + E[3] * sqrt(ρᵢ)
    Fl = F[1] + F[2] * log(ρᵢ) + F[3] * exp(-ρᵢ)
    Gl = (G[1] + G[2] * log(ρᵢ) * sqrt(ρᵢ) + G[3] / sqrt(ρᵢ))^(-1)
    Hl = H[1] + H[2] * (ρᵢ)^(5 / 2) - H[3] * exp(-ρᵢ)
    # Table B4
    ai = (Bl * ρ^Al, El * ρ^Al * exp(Hl * ρ))
    bi = (Cl, Fl)
    ci = (FT(0), Gl)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return (aiu, bi, ciu)
end

# Wrapper to cast types from SF.gamma
# (which returns Float64 even when the input is Float32)
# TODO - replace with parameterization of our own
Γ(a::FT) where {FT <: Real} = FT(SF.gamma(a))

"""
    Chen2022_monodisperse_pdf(a, b, c, D)

 - a, b, c, - free parameters defined in Chen etl al 2022
 - D - droplet diameter

Returns the addends of the bulk fall speed of rain or ice particles
following Chen et al 2022 DOI: 10.1016/j.atmosres.2022.106171 in [m/s].
Assuming monodisperse droplet distribution.
"""
function Chen2022_monodisperse_pdf(a, b, c, D)
    return a * D^b * exp(-c * D)
end

"""
    Chen2022_exponential_pdf(a, b, c, λ, k)

 - a, b, c, - free parameters defined in Chen etl al 2022
 - λ - size distribution parameter
 - k - size distribution moment for which we compute the bulk fall speed

Returns the addends of the bulk fall speed of rain or ice particles
following Chen et al 2022 DOI: 10.1016/j.atmosres.2022.106171 in [m/s].
Assuming exponential size distribution and hence μ=0.
"""
function Chen2022_exponential_pdf(a::FT, b::FT, c::FT, λ::FT, k::Int) where {FT}
    μ = 0 # Exponential instead of gamma distribution
    δ = FT(μ + k + 1)
    return a * λ^δ * Γ(b + δ) / (λ + c)^(b + δ) / Γ(δ)
end

end # module end

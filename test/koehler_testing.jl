using Plots

# Atmospheric conditions
T = 294.0         # air temperature
p = 1000.0 *1e2   # air pressure
w = 0.5           # vertical velocity

# NaCl - universal parameters
M_s = 0.05844     # molar mass
ρ_s = 2.16 * 1e3  # density
ϕ_s = 1.0         # osmotic coeff
ν_s = 2.0         # n ions

# kappa
κ = 1.28

# water
M_l = 0.01801     # molar mass
ρ_l = 1e3         # density
σ = 0.072         # surface tension

# consts
R = 8.3145
Rv = R / M_l

# dry radius
rd = 0.02 * 1e-6

function A_fun(T)
    2 * σ / ρ_l / Rv / T
end

function B_fun(rd)
    ν_s * ϕ_s * rd^3 * ρ_s * M_l / ρ_l / M_s
end

function S_1(r, rd, T)
    A = A_fun(T)
    B = B_fun(rd)
    exp.(A./r - B./r.^3)
end

function S_1_approx(r, rd, T)
    A = A_fun(T)
    B = B_fun(rd)
    1.0 .+ A./r - B./r.^3
end

function S_1_A_term(r, T)
    A = A_fun(T)
    1.0 .+ A./r
end

function S_1_B_term(r, rd)
    B = B_fun(rd)
    1.0 .- B./r.^3
end

function aw_fun(r, rd)
    Vs = 4.0/3 * 3.14 * rd^3
    Vl = 4.0/3 * 3.14 .*(r.^3 .- rd^3)
    1.0 ./ (1.0 .+ κ .* Vs ./ Vl)
end

function S_2(r, rd, T)
    aw = aw_fun(r, rd)
    A = A_fun(T)
    aw .* exp.(A ./ r)
end

function S_2_approx(r, rd, T)
    A = A_fun(T)
    (r.^3 .- rd^3) ./ (r.^3 .- rd^3*(1.0 - κ)) .* exp.(A ./ r)
end

r_range = range(0.1 * 1e-6, stop=10 * 1e-6, length=1000)

#plot(r_range * 1e6, S_1(r_range, rd, T), label="S_1", xaxis=(:log10), xlabel="rw [um]", ylabel="S")
#plot!(r_range * 1e6, S_1_approx(r_range, rd, T), label="S_1_approx", xaxis=(:log10), xlabel="rw [um]", ylabel="S")

#plot!(r_range * 1e6, S_1_A_term(r_range, T), label="S_1_A_term", xaxis=(:log10), xlabel="rw [um]", ylabel="S")
#plot!(r_range * 1e6, S_1_B_term(r_range, rd), label="S_1_B_term", xaxis=(:log10), xlabel="rw [um]", ylabel="S")

plot(r_range * 1e6, S_2(r_range, rd, T), label="S_2", xaxis=(:log10), xlabel="rw [um]", ylabel="S")
plot!(r_range * 1e6, S_2(r_range, rd, T), label="S_2_approx", xaxis=(:log10), xlabel="rw [um]", ylabel="S")

savefig("Koehler.svg")

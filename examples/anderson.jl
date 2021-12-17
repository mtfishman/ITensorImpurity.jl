using ITensors
using ITensorImpurity
using UnicodePlots
using Random

Random.seed!(1234)

# Number of sites
# N = 16, 32, 64, 128
N = 8
Nᴸ = N ÷ 2
Nᴿ = N ÷ 2 - 1
@assert iseven(Nᴸ)
@assert isodd(Nᴿ)
@assert Nᴸ + Nᴿ + 1 == N

@show N
@show Nᴸ
@show Nᴿ

# Bath hopping
t = 1.0

# Time steps
# Typical time steps: δτ = 0.01 - 0.1
# Can use larger time steps, up to δτ ≈ 0.4,
# deeper in the Kondo regime (large U/t′²)
# Reference: https://arxiv.org/pdf/0807.0581.pdf
δτ = 0.05t
τf = 6.0t

# Hybridization paramater
U = t
t′ = 0.4
# Source: https://arxiv.org/pdf/0807.0581.pdf
# U / Γ = 3.125
Γ = 2 * t′ ^ 2 / t

# Bias potential
V = 0.005U

# Impurity potential
# Particle-hole symmetric point is: Vᵍ = -U / 2
Vᵍ = -U / 2

@show t, U, Γ, U / Γ, t′, V, Vᵍ

s = siteinds("Electron", N)
# Without bias potential
ℋ₀ = anderson(N; Nᴸ, t, V=0.0, t′, U, Vᵍ)
# With bias potential
ℋ = anderson(N; Nᴸ, t, V, t′, U, Vᵍ)

H₀ = MPO(ℋ₀, s)
H = MPO(ℋ, s)

ψ0 = randomMPS(s; linkdims=4)

sweeps = Sweeps(10)
setmaxdim!(sweeps, 50, 100)
setcutoff!(sweeps, 1e-8)
e, ψ0 = dmrg(H₀, ψ0, sweeps)

ntot = expect(ψ0, "ntot")
println("⟨ψ|ntot|ψ⟩ = ")
display(ntot)
@show N
@show sum(ntot)

expH = exp(Val{:trotter_order_2}(), -im * δτ, ℋ, s)
U = MPO(s, "Id")
U = apply(expH, U; cutoff=1e-6)

if N ≤ 4
  # Compare trotter decomposition with ED
  ∏U = contract_slice(linkinds, U)
  ∏H = contract_slice(linkinds, H)
  @show norm(∏U - exp(-im * δt * ∏H))
end

τi = 0.0t
τ⃗ = τi:δτ:τf
nlayers = length(τ⃗)
f(ψτ) = ITensorImpurity.J(ψτ, t′; Nᴸ) / V
ψτ, J⃗ = apply_layers([expH for _ in 1:(nlayers - 1)], ψ0; f, cutoff=1e-8)

p = lineplot(τ⃗, real(J⃗))
display(p)

nothing

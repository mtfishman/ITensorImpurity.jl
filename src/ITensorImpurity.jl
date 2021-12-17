module ITensorImpurity
  using ITensors

  export anderson, apply_layers

  function ITensors.op(o::ITensors.SiteOp, s::AbstractVector{<:Index})
    return op(o.name, s, o.site)
  end

  ITensors.sites(𝒽::ITensors.MPOTerm) = [only(o.site) for o in 𝒽.ops]
  ITensors.names(𝒽::ITensors.MPOTerm) = [o.name for o in 𝒽.ops]

  function ITensors.op(𝒽::ITensors.MPOTerm, s::AbstractVector{<:Index})
    n⃗ = ITensors.sites(𝒽)
    p = sortperm(n⃗)
    opnames = ITensors.names(𝒽)
    s⃗ = s[n⃗][p]
    𝒽′ = ITensors.MPOTerm(𝒽.coef, [ITensors.SiteOp(opnames[n], p[n]) for n in eachindex(𝒽.ops)])
    return prod(MPO(ITensors.OpSum([𝒽′]), s⃗))
  end

  # Trotter layer
  function Base.exp(::Val{:trotter_order_1}, δt::Number, ℋ::OpSum, s::AbstractVector{<:Index})
    # TODO: Sort and merge terms on the same sites.
    circuit = Vector{ITensor}(undef, length(ℋ))
    for n in 1:length(ℋ)
      # TODO: Use `cis` in the imaginary `δt` case.
      circuit[n] = exp(δt * op(ℋ[n], s))
    end
    return circuit
  end

  function Base.exp(::Val{:trotter_order_2}, δt::Number, ℋ::OpSum, s::Vector{<:Index})
    expH1 = exp(Val{:trotter_order_1}(), δt / 2, ℋ, s)
    return vcat(expH1, reverse(expH1))
  end

  function contract_slice(t1::ITensor, t2::ITensor; slice_ind)
    @assert hasind(t1, slice_ind)
    @assert hasind(t2, slice_ind)
    i = slice_ind
    r = ITensor()
    for n in 1:dim(i)
      pn = onehot(i => n)
      r += (t1 * pn) * (t2 * pn)
    end
    return r
  end

  function contract_slice(::typeof(linkinds), tn)
    ∏tn = tn[1]
    for n in 2:length(tn)
      l = commonind(∏tn, tn[n])
      ∏tn = contract_slice(∏tn, tn[n]; slice_ind=l)
    end
    return ∏tn
  end

  const σ⃗ = ("↑", "↓")

  function anderson(N; Nᴸ, t, V, t′, U, Vᵍ)
    ℋ = OpSum()
    for n in 1:Nᴸ
      if n < Nᴸ
        for σ in σ⃗
          ℋ .-= t, "c†$σ", n, "c$σ", n + 1
          ℋ .-= t, "c†$σ", n + 1, "c$σ", n
        end
      end
      ℋ .-= V / 2, "ntot", n
    end
    for σ in σ⃗
      ℋ .-= t′, "c†$σ", Nᴸ, "c$σ", Nᴸ + 1
      ℋ .-= t′, "c†$σ", Nᴸ + 1, "c$σ", Nᴸ
      ℋ .-= t′, "c†$σ", Nᴸ + 1, "c$σ", Nᴸ + 2
      ℋ .-= t′, "c†$σ", Nᴸ + 2, "c$σ", Nᴸ + 1
    end
    ℋ .+= U, "n↑↓", Nᴸ + 1
    ℋ .+= Vᵍ, "ntot", Nᴸ + 1
    for n in (Nᴸ + 2):N
      if n < N
        for σ in σ⃗
          ℋ .-= t, "c†$σ", n, "c$σ", n + 1
          ℋ .-= t, "c†$σ", n + 1, "c$σ", n
        end
      end
      ℋ .+= V / 2, "ntot", n
    end
    return ℋ
  end

  function J(ψᵗ, t′; Nᴸ)
    J⃗ᴸᵗ = Dict()
    J⃗ᴿᵗ = Dict()
    for σ in σ⃗
      J⃗ᴸᵗ[σ] = im * t′ * correlation_matrix(ψᵗ, "c†$σ", "c$σ"; site_range=Nᴸ:(Nᴸ + 1))
      J⃗ᴿᵗ[σ] = im * t′ * correlation_matrix(ψᵗ, "c†$σ", "c$σ"; site_range=(Nᴸ + 1):(Nᴸ + 2))
    end
    Jᴸᵗ = sum(σ -> J⃗ᴸᵗ[σ][1, 2] - J⃗ᴸᵗ[σ][2, 1], σ⃗)
    Jᴿᵗ = sum(σ -> J⃗ᴿᵗ[σ][1, 2] - J⃗ᴿᵗ[σ][2, 1], σ⃗)
    return (Jᴸᵗ + Jᴿᵗ) / 2
  end

  function cc(ψ; Nᴸ, t′)
    ccs = Dict()
    for σ in σ⃗
      ccs[σ] = im * t′ * correlation_matrix(ψ, "c†$σ", "c$σ"; site_range=Nᴸ:(Nᴸ + 2))
    end
    return ccs
  end

  function apply_layers(U⃗, ψ0; f, kwargs...)
    nlayers = length(U⃗)
    f⃗ = Vector{Any}(undef, nlayers + 1)
    ψτ = ψ0
    f⃗[1] = f(ψτ)
    for n in 1:nlayers
      @show n, nlayers
      ψτ = apply(U⃗[n], ψτ; kwargs...)
      @show maxlinkdim(ψτ)
      @show fτ = f(ψτ)
      f⃗[n + 1] = fτ
    end
    return ψτ, f⃗
  end
end # module ITensorImpurity

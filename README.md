# ITensorImpurity

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mtfishman.github.io/ITensorImpurity.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mtfishman.github.io/ITensorImpurity.jl/dev)
[![Build Status](https://github.com/mtfishman/ITensorImpurity.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mtfishman/ITensorImpurity.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mtfishman/ITensorImpurity.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mtfishman/ITensorImpurity.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Installation

This package is not yet registered, install it as follows:
```julia
julia> using Pkg

julia> Pkg.add("https://github.com/mtfishman/ITensorImpurity.jl")
```
If you want the development version to make local changes, you can then do:
```julia
julia> Pkg.develop("ITensorImpurity")
```
This will make a clone of this repository in the directory `~/.julia/dev/ITensorImpurity`, which (in conjunction with [Revise.jl](https://github.com/timholy/Revise.jl)) will allow you to make local changes to the package and modify the example code.

## Overview

Currently, this is a small expiremental package that provides a definition of an Anderson impurity mode Hamiltonian and some tools for time evolving it.

Once installed, the main interface is the function `anderson` which outputs an Andersion impurity Hamiltonian:
```julia
# Number of sites
N = 8

# Location of the impurity
Nᴸ = N ÷ 2

# Bath hopping
t = 1.0

# Hybridization paramater
U = t

# Impurity hopping
t′ = 0.4

# Bias potential
V = 0.005U

# Impurity potential
# Particle-hole symmetric point is: Vᵍ = -U / 2
Vᵍ = -U / 2

# Representation of the Hamiltonian
ℋ = anderson(N; Nᴸ, t, V, t′, U, Vᵍ)

# Convert the Hamiltonian to an MPO
s = siteinds("Electron", N)
H = MPO(ℋ, s)
```
For a full example of running DMRG and time evolving with a bias potential, see the [examples](https://github.com/mtfishman/ITensorImpurity.jl/tree/main/examples).

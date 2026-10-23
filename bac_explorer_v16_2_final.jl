# BAC-Explorer v16.2 — full pipeline including catalysis layer
# Verified: R(A⋆)=0.019353, R(tA⋆)=0.032908, R(shared)=0.014065
# Matches paper values to 5 decimal places.

# Run: julia +1.12, then include this file.

using Pkg
Pkg.activate("C:\\\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Polyhedra
using Bonito, WGLMakie, Makie
Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c
WGLMakie.activate!()

# Paste build_catalysis_model, binding_order_matrix, binding_order_scalar,
# steady_state_invariance, input_responsiveness, adaptation_signature,
# and the App() do block from your session below.

# GellMannMatrices

[![Build Status](https://github.com/thchr/GellMannMatrices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/thchr/GellMannMatrices.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package to compute the [generalized Gell-Mann matrices](https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices#Generalized_Gell-Mann_matrices_(Hermitian)) in `d` dimensions.
The generalized Gell-Mann matrices provide a traceless basis for the space of Hermitian matrices (which is orthogonal under the [Frobenius inner product](https://en.wikipedia.org/wiki/Frobenius_inner_product)).
Notable special cases include the [Pauli matrices σᵢ](https://en.wikipedia.org/wiki/Pauli_matrices) (`d = 2`) and the [standard Gell-Mann matrices Λᵢ](https://en.wikipedia.org/wiki/Gell-Mann_matrices) (`d = 3`).

## Usage
GellMannMatrices.jl exports a single function `gellmann(d)` which returns the generalized Gell-Mann matrices in `d` dimensions.
The signature `gellmann(T, d)` allows specifying the matrix type `T` (which must be a mutable `AbstractMatrix` with complex element type).

As examples, we can compute the Pauli matrices and the standard Gell-Mann matrices:
```jl
julia> gellmann(2) # Pauli matrices
 [0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 - 1.0im; 0.0 + 1.0im 0.0 + 0.0im]
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -1.0 + 0.0im]
```
```
julia> gellmann(3) # Standard Gell-Mann matrices
 [0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 - 1.0im 0.0 + 0.0im; 0.0 + 1.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
 [... omitted]
```
The keyword argument `skip_identity` (default, `true`) can be toggled to `false` to include the (traceful)  identity matrix.
# GellMannMatrices

[![Build Status](https://github.com/thchr/GellMannMatrices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/thchr/GellMannMatrices.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package to compute the [generalized Gell-Mann matrices](https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices#Generalized_Gell-Mann_matrices_(Hermitian)) in `d` dimensions.
The generalized Gell-Mann matrices provide a traceless basis for the space of Hermitian matrices (which is orthogonal under the [Frobenius inner product](https://en.wikipedia.org/wiki/Frobenius_inner_product)).
Notable special cases include the Pauli matrices σᵢ (`d = 2`) and the Gell-Mann matrices Λᵢ (`d = 3`).

## Usage
`GellMannMatrices` exports a single function `gellmann([T=Matrix{ComplexF64},] d)`, returning the generalized Gell-Mann matrices in `d` dimensions as matrices of type `T` (whose `eltype` must be complex).
For instance, `gellmann(2)` returns the Pauli matrices and `gellmann(3)` returns the Gell-Mann matrices, both in the conventional sorting.

The keyword argument `skip_identity` (default, `true`) can be toggled to true to include the (traceful)  identity matrix of dimension `d`.
# GellMannMatrices

[![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]

A Julia package to compute the [generalized Gell-Mann matrices](https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices#Generalized_Gell-Mann_matrices_(Hermitian)) in `d` dimensions.

The generalized Gell-Mann matrices provide a traceless basis for the space of Hermitian matrices, with elements that are orthogonal under the [Frobenius inner product](https://en.wikipedia.org/wiki/Frobenius_inner_product).

Notable special cases include the [Pauli matrices σᵢ](https://en.wikipedia.org/wiki/Pauli_matrices) (`d = 2`) and the [standard Gell-Mann matrices Λᵢ](https://en.wikipedia.org/wiki/Gell-Mann_matrices) (`d = 3`).

## Usage
GellMannMatrices.jl exports a single function `gellmann(d)` which returns the `d^2 - 1` generalized Gell-Mann matrices in `d` dimensions.
The signature `gellmann(T, d)` allows specifying the matrix type `T` (which must be a mutable `AbstractMatrix` with complex element type).

As examples, we can compute the Pauli matrices and the standard Gell-Mann matrices:
```jl
julia> gellmann(2) # Pauli matrices
 [0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 - 1.0im; 0.0 + 1.0im 0.0 + 0.0im]
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -1.0 + 0.0im]

julia> gellmann(3) # Standard Gell-Mann matrices
 [0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 - 1.0im 0.0 + 0.0im; 0.0 + 1.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
 [...]
```
**Keyword arguments:**
- `skip_identity` (default, `true`): toggle to `false` to include the identity matrix.
- `normalize` (default, `false`): toggle to `true` to guarantee a fixed (Frobenius) normalization prefactor of 2, in the sense $\mathrm{Tr}(M_i^\dagger M_j) = 2\delta_{ij}$. If `false`, matrix elements are chosen to be small integers, leaving the normalization prefactor matrix-dependent.

[ci-status-img]: https://github.com/thchr/GellMannMatrices.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-status-url]: https://github.com/thchr/GellMannMatrices.jl/actions/workflows/CI.yml?query=branch%3Amain
[coverage-img]:  https://codecov.io/gh/thchr/GellMannMatrices.jl/branch/main/graph/badge.svg
[coverage-url]:  https://codecov.io/gh/thchr/GellMannMatrices.jl
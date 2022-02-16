module GellMannMatrices

using LinearAlgebra, SparseArrays
export gellmann

# ---------------------------------------------------------------------------------------- #
# until we have https://github.com/JuliaLang/julia/pull/44211, we need to extend `dims2cat`
# to ensure type-stability of `cat(..., dims=(1,2))`
# ... we probably shouldn't do this, piracy and all that :(
function Base.dims2cat(::Val{dims}) where dims
    if any(≤(0), dims)
        throw(ArgumentError("All cat dimensions must be positive integers, but got $dims"))
    end
    ntuple(in(dims), maximum(dims))
end

# ---------------------------------------------------------------------------------------- #

# direct sum (i.e., block-diagonalization of matrices As)
@inline ⊕(As...) = cat(As..., dims=Val((1,2)))

# https://en.wikipedia.org/wiki/Matrix_unit
unit(::Type{<:SparseMatrixCSC}, i, j, n, m) = sparse([i], [j], [1], n, m)
unit(T::Type{<:AbstractMatrix}, i, j, n, m) = (A = zeros(eltype(T), n, m); A[i,j] = one(eltype(T)); A)

# ---------------------------------------------------------------------------------------- #
# implementation following https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices

# off-diagonal matrices (f_ij on Wiki)
function offdiag_gellmann(T::Type{<:AbstractMatrix}, i, j, d)
    if i < j
        return unit(T, i, j, d, d) + unit(T, j, i, d, d)
    elseif i > j
        return -1im * (unit(T, j, i, d, d) - unit(T, i, j, d, d))
    else
        throw(DomainError((j, i), "i and j must differ"))
    end
end

# "diagonal"-like matrices (h_i on Wiki)
function diag_gellmann(T::Type{<:AbstractMatrix}, i, d)
    if d == 0
        return T <: SparseMatrixCSC ? spzeros(eltype(T), 0, 0) :
               T <: AbstractMatrix ?  zeros(eltype(T), 0, 0)   : error()
    end
    if i == 1
        T(I(d))
    elseif i < d
        diag_gellmann(T, i, d-1) ⊕ zero(eltype(T))
    elseif i == d
        sqrt(2/(d*(d-1))) * (T(I(d-1)) ⊕ convert(eltype(T), (1-d)))
    else
        error()
    end
end

# ---------------------------------------------------------------------------------------- #

"""
    gellmann([T::Type{<:AbstractMatrix}=Matrix{ComplexF64},] d::Integer; skip_identity=true)

Return the generalized Gell-Mann matrices in `d` dimensions as a vector of matrices of type
`T`.

`T` can be any mutable `AbstractMatrix`, including `SparseMatrixCSC`, with complex floating
point element type.

## Keyword arguments
- `skip_identity` (default, `false`): if `true`, the identity matrix of dimension `d` is
included as the last element of the returned vector.

## Examples
Compute the Pauli matrices ``σ_i`` (`d=2`) and Gell-Mann matrices ``Λ_i`` (`d=3`):
```jl
julia> σᵢ = gellmann(2)
julia> Λᵢ = gellmann(3)
```
"""
gellmann(d::Integer; kwargs...) = gellmann(Matrix{ComplexF64}, d; kwargs...)

function gellmann(T::Type{<:AbstractMatrix}, d::Integer; skip_identity=true)
    eltype(T) == complex(eltype(T)) || throw(DomainError(T, "the element type of T must be complex"))
    # we collect the matrices in the slightly odd-looking manner below in order to be
    # consistent with the conventional number of Pauli and Gell-Mann matrices.
    ms = Vector{T}(undef, d*d - skip_identity)
    k = 0
    for i in 2:d
        for  j in 1:i-1
            ms[k+=1] = offdiag_gellmann(T, j, i, d)
            ms[k+=1] = offdiag_gellmann(T, i, j, d)
        end
        ms[k+=1] = diag_gellmann(T, i, d)
    end
    skip_identity || (ms[k+=1] = diag_gellmann(T, 1, d))
    return ms
end

end

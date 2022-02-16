using GellMannMatrices
using LinearAlgebra
using SparseArrays
using Test

@testset "GellMannMatrices.jl" begin
    # d=2
    @test (@inferred gellmann(1)) == []
    @test (@inferred gellmann(SparseMatrixCSC{ComplexF64, Int}, 1)) == []
    @test (@inferred gellmann(1; skip_identity=false)) == [ones(1,1)]

    # Pauli matrices (d=2)
    σᵢ = [[0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1]]
    @test (@inferred gellmann(2)) ≈ σᵢ
    @test (@inferred gellmann(SparseMatrixCSC{ComplexF64, Int}, 2)) ≈ σᵢ
    @test (@inferred gellmann(2; skip_identity=false)) ≈ [σᵢ..., [1 0; 0 1]]

    # Gell-Mann matrices (d=3)
    Λᵢ = [[0 1 0; 1 0 0; 0 0 0], [0 -1im 0; 1im 0 0; 0 0 0], [1 0 0; 0 -1 0; 0 0 0],
          [0 0 1; 0 0 0; 1 0 0], [0 0 -1im; 0 0 0; 1im 0 0], [0 0 0; 0 0 1; 0 1 0],
          [0 0 0; 0 0 -1im; 0 1im 0], [1 0 0; 0 1 0; 0 0 -2]/√3]
    @test (@inferred gellmann(3)) ≈ Λᵢ
    @test (@inferred gellmann(SparseMatrixCSC{ComplexF64, Int}, 3)) ≈ Λᵢ
    @test (@inferred gellmann(3; skip_identity=false)) ≈ [Λᵢ..., [1 0 0; 0 1 0; 0 0 1]]

    # hermicity, tracelessness, & orthogonality
    for d in 1:3
        ms = gellmann(d)
        for (i,mᵢ) in enumerate(ms)
            @test mᵢ ≈ mᵢ'
            @test norm(tr(mᵢ)) < 1e-12
            for (j,mⱼ) in enumerate(ms)
                if i == j
                    @test abs(tr(mᵢ'*mᵢ)) ≈ 2 # Frobenius matrix norm
                else
                    @test abs(tr(mᵢ'*mⱼ)) < 1e-12
                end
            end
        end
    end
end

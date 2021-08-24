# Test Gaussian tetrahedron integration
using TetrahedronIntegration
using Test

using LinearAlgebra
using HCubature
using StaticArrays
using TetrahedronIntegration: integrate_gaussian_times_polynomial

@testset "gaussian" begin
    @testset "gaussian times polynomial" begin
        σ = 0.7
        coeffs = (0.3, 0.4, 0.7, 1.1)
        a = -0.9
        b = 1.2
        f(x) = exp(-x^2/σ^2) * evalpoly(x, coeffs)

        @test integrate_gaussian_times_polynomial(σ, coeffs, a, b) ≈ hquadrature(f, a, b)[1]
    end

    @testset "gaussian parallelepiped" begin
        f(x) = exp(-(e0 + dot(x, v0))^2 / σ^2) / sqrt(π)
        σ = 0.7
        e0 = 0.3
        v0 = SVector(0.5, -0.3, 0.2)
        L = 0.4
        Lvec =  SVector((L, L, L))
        @test gaussian_parallelepiped(σ, e0, v0, L) ≈ hcubature(f, -Lvec/2, Lvec/2)[1] ./ L^3

        v0 = SVector(0.0, 0.0, 0.1)
        @test gaussian_parallelepiped(σ, e0, v0, L) ≈ hcubature(f, -Lvec/2, Lvec/2)[1] ./ L^3

        v0 = SVector(0.0, 0.0, 0.0)
        @test gaussian_parallelepiped(σ, e0, v0, L) ≈ hcubature(f, -Lvec/2, Lvec/2)[1] ./ L^3
    end
end
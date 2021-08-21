# Test tetrahedron delta function integration
using TetrahedronIntegration
using StaticArrays
using Test
using TetrahedronIntegration: _evaluate_piecewise_polynomial

@testset "delta" begin
    e0 = 0.4
    v0 = [0.0, 0.0, 1.0]
    L = 0.3
    de = L * v0[3] / 2
    @test delta_parallelepiped(e0, e0, v0, L) ≈ 1.0 / L
    @test delta_parallelepiped(e0 + de - 1e-3, e0, v0, L) ≈ 1.0 / L
    @test delta_parallelepiped(e0 - de + 1e-3, e0, v0, L) ≈ 1.0 / L
    @test delta_parallelepiped(e0 + de + 1e-3, e0, v0, L) ≈ 0
    @test delta_parallelepiped(e0 - de - 1e-3, e0, v0, L) ≈ 0

    # Test small v0
    @test delta_parallelepiped(e0, e0, [0.0, 0.0, 0.0], L) ≈ 0
end

@testset "delta_polynomial" begin
    e1234 = SVector(1.5, -0.2, 0.7, 0.3)
    polys = delta_tetrahedron_polynomial(e1234)
    for e in -0.3:0.1:1.6
        @test delta_tetrahedron(e, e1234) ≈ _evaluate_piecewise_polynomial(e, polys) atol=10*eps(eltype(e))
    end

    e1234 = SVector(1.0, -0.5, 1.0, -0.5)
    polys = delta_tetrahedron_polynomial(e1234)
    for e in -1.0:0.5:1.0
        @test delta_tetrahedron(e, e1234) ≈ _evaluate_piecewise_polynomial(e, polys)
    end
end
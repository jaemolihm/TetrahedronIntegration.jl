using LinearAlgebra
using Polynomials
using StaticArrays
using SpecialFunctions

export gaussian_parallelepiped

"""
Compute ∫_a^b dx exp(-x^2/σ^2) * poly(x) = σ * ∫_{a/σ}^{b/σ} dx exp(-x^2) * poly(σx)
"""
function integrate_gaussian_times_polynomial(σ::FT, poly::AbstractPolynomial, a, b) where {FT}
    length(poly.coeffs) > 4 && error("Quartic polynomial not implemented")
    aa = a / σ
    bb = b / σ
    sqrt_pi = sqrt(typeof(σ)(π))
    erf_ba = erf(aa, bb)
    exp2_a = exp(-aa^2)
    exp2_b = exp(-bb^2)
    val = FT(0)
    for (i, c) in enumerate(poly.coeffs)
        if i == 1
            val += c * sqrt_pi * erf_ba / 2
        elseif i == 2
            val += -c * σ * (exp2_b - exp2_a) / 2
        elseif i == 3
            val += c * σ^2 * (2aa * exp2_a - 2bb * exp2_b + sqrt_pi * erf_ba) / 4
        elseif i == 4
            val += c * σ^3 * ((aa^2 + 1) * exp2_a - (bb^2 + 1) * exp2_b) / 2
        end
    end
    σ * val
end

"""
Calculate ``1/volume * ∫_{tetrahedron} d^3x exp(-e(x)^2/σ^2)`` where e1234 are the
values of e(x) at the four vertices of the tetrahedron.
"""
@inline function gaussian_tetrahedron(σ::FT, e1234) where {FT}
    polys = delta_tetrahedron_polynomial(e1234)
    val = zero(FT)
    for (e1, e2, poly) in polys
        e1 == e2 && continue
        val += integrate_gaussian_times_polynomial(σ, poly, e1, e2)
    end
    val
end

"""
Calculate
``1/volume * ∫_{parallelepiped} d^3x f(e(x))
= 1/volume * ∫_{parallelepiped} ∫ dee d^3x f(ee) delta(ee - e(x))
= ∫ dee f(ee) tetra(ee)``,
where ``f(ee) = exp(-ee^2/σ^2)``,
``tetra(ee) = 1/volume * ∫_{parallelepiped} d^3x delta(ee - e(x))``,
and e(x) = `e0` + x ⋅ `v0`. The parallelepiped region is [-L/2, L/2]^3.
Divide the parallelepiped into six tetrahedra and use tetrahedron integration.
"""
function gaussian_parallelepiped(σ::FT, e0, v0, L) where {FT}
    if L * norm(v0) < σ * 1e-4
        # Velocity is too small. Tetrahedron is unstable, use single point.
        return exp(-e0^2/σ^2)
    end
    L_div_2 = L ./ 2
    v0_L = v0 .* L_div_2
    e1 = e0 - v0_L[1] - v0_L[2] - v0_L[3]
    e2 = e0 + v0_L[1] - v0_L[2] - v0_L[3]
    e3 = e0 - v0_L[1] + v0_L[2] - v0_L[3]
    e4 = e0 + v0_L[1] + v0_L[2] - v0_L[3]
    e5 = e0 - v0_L[1] - v0_L[2] + v0_L[3]
    e6 = e0 + v0_L[1] - v0_L[2] + v0_L[3]
    e7 = e0 - v0_L[1] + v0_L[2] + v0_L[3]
    e8 = e0 + v0_L[1] + v0_L[2] + v0_L[3]
    val = zero(FT)
    val += gaussian_tetrahedron(σ, SVector{4,FT}(e1, e8, e6, e2))
    val += gaussian_tetrahedron(σ, SVector{4,FT}(e1, e8, e2, e4))
    val += gaussian_tetrahedron(σ, SVector{4,FT}(e1, e8, e4, e3))
    val += gaussian_tetrahedron(σ, SVector{4,FT}(e1, e8, e3, e7))
    val += gaussian_tetrahedron(σ, SVector{4,FT}(e1, e8, e7, e5))
    val += gaussian_tetrahedron(σ, SVector{4,FT}(e1, e8, e5, e6))
    return val / 6
end
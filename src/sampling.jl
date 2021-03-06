using StaticArrays

export delta_parallelepiped_sampling!

"""
Calcultate val = int_{parallelepiped} d^3k delta(`etarget` - e(k)) / volume
where `e(k) = e0 + k ⋅ v0`. The parallelepiped region is [-L/2, L/2]^3.
Divide the parallelepiped into six tetrahedra and use tetrahedron integration.
# Sampling
Sample points on the `e(k) = etarget` plane, where the spacing is relative to the edge
lengths. Sample around `nsample_1d` points along each dimension, and at most
`(nsample_1d+1)^2` points are sampled. Sampled points are stored in the array `ksamples`.
"""
@inline function delta_parallelepiped_sampling!(etarget::T, e0, v0, L, nsample_1d, ksamples) where {T}
    if norm(v0) < 1E-10
        # Velocity is too small. Tetrahedron method does not work.
        return zero(T), 0
    end

    L_div_2 = L ./ 2
    e1 = e0 + (-v0[1] * L_div_2[1] - v0[2] * L_div_2[2] - v0[3] * L_div_2[3])
    e2 = e0 + (+v0[1] * L_div_2[1] - v0[2] * L_div_2[2] - v0[3] * L_div_2[3])
    e3 = e0 + (-v0[1] * L_div_2[1] + v0[2] * L_div_2[2] - v0[3] * L_div_2[3])
    e4 = e0 + (+v0[1] * L_div_2[1] + v0[2] * L_div_2[2] - v0[3] * L_div_2[3])
    e5 = e0 + (-v0[1] * L_div_2[1] - v0[2] * L_div_2[2] + v0[3] * L_div_2[3])
    e6 = e0 + (+v0[1] * L_div_2[1] - v0[2] * L_div_2[2] + v0[3] * L_div_2[3])
    e7 = e0 + (-v0[1] * L_div_2[1] + v0[2] * L_div_2[2] + v0[3] * L_div_2[3])
    e8 = e0 + (+v0[1] * L_div_2[1] + v0[2] * L_div_2[2] + v0[3] * L_div_2[3])
    val = zero(T)
    val += delta_tetrahedron(etarget, SVector{4,T}(e1, e2, e4, e5))
    val += delta_tetrahedron(etarget, SVector{4,T}(e1, e3, e4, e5))
    val += delta_tetrahedron(etarget, SVector{4,T}(e2, e4, e6, e5))
    val += delta_tetrahedron(etarget, SVector{4,T}(e3, e4, e7, e5))
    val += delta_tetrahedron(etarget, SVector{4,T}(e4, e7, e8, e5))
    val += delta_tetrahedron(etarget, SVector{4,T}(e4, e6, e8, e5))

    # TODO: If |v0| = 0, skip sampling

    # Find one point on the e(k0) = etarget plane
    e = (e1, e2, e3, e4, e5, e6, e7, e8)
    imin = argmin(e)
    imax = argmax(e)
    emin = e[imin]
    emax = e[imax]
    # Trivial case: etarget outside the energy range
    if etarget < emin || emax < etarget
        return zero(T), 0
    end

    kmin, kmax = KLIST[imin], KLIST[imax]
    ratio = (etarget - emin) / (emax - emin)
    k0 = (kmin + (kmax - kmin) * ratio) .* L_div_2

    # Find two orthogonal unit vectors on the plane
    i = argmin(abs.(v0))
    if i == 1
        u0 = Vec3{T}(1, 0, 0)
    elseif i == 2
        u0 = Vec3{T}(0, 1, 0)
    else
        u0 = Vec3{T}(0, 0, 1)
    end
    u1 = cross(v0, u0)
    u2 = cross(v0, u1)

    # Normalize unit vectors according to the step size
    u1 = u1 / norm(u1 ./ L) / nsample_1d
    u2 = u2 / norm(u2 ./ L) / nsample_1d

    # k vectors on the e(k) = etarget plane is k0 + n1 * u1 + n2 * u2.
    # Sample k vectors on the plane and inside the cube.

    nmax = floor(Int, sqrt(3) * nsample_1d)
    nsamples = 0
    for n1 in -nmax:nmax
        was_inside = false
        k = k0 + n1 * u1 - nmax * u2
        for n2 in -nmax:nmax
            k += u2
            if all(abs.(k) .<= L_div_2)
                was_inside = true
                nsamples += 1
                ksamples[:, nsamples] .= k
            elseif was_inside
                # Was inside the cube, but moved out. The remaining points are all outside.
                break
            end
        end
    end

    return val / 6, nsamples
end
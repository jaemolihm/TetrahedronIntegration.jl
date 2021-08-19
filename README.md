# TetrahedronIntegration.jl

[![Build Status](https://github.com/jaemolihm/TetrahedronIntegration.jl/workflows/CI/badge.svg)](https://github.com/jaemolihm/TetrahedronIntegration.jl/actions)
[![Coverage](https://codecov.io/gh/jaemolihm/TetrahedronIntegration.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jaemolihm/TetrahedronIntegration.jl)

Tetrahedron integration


## References
* P. B. Allen, Phys. Stat. Sol. (b) 120, 629 (1983)
* Blöchl et al, Phys. Rev. B 16, 232 (1994)

## TODO
* Add introduction to README
* Add references
* Full coverage
* Docs
* Improve gaussian with small e2-e1, e3-e2, e4-e3, Add tests
* Implement Blöchl correction (Eq. (22))
* Which diagonal to cut the cube into tetrahedron
* Modularize (e0, v0) -> (e1, ..., e8)
* Add type for piecewise polynomial (?)
* Functions other than gaussian (by quadrature integration?)
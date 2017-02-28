Conic programming library for quantum information
=================================================

This library provides several useful cones for quantum information research.

- `SeparableConeC` provides exact/inner and outer approximations of the cone
  of bipartite separable operators

- `AbsoluteRobustnessConeC`, `GeneralizedRobustnessConeC`, `RandomRobustnessConeC`
  provide formulations of various robustness-based entanglement measures

- `DUBConeC` provides an entanglement cone formulation of an entanglement measure
  due to Duan&Wang, 2016

- `NegativityConeC` provides a conic formulation of the negativity entanglement measure

Some classes/functions are used internally by the library, but will be useful in
other contexts:

- `SymmetricSubspace` provides methods to work with a canonical basis of the symmetric
  subspace: index <-> subindices conversion, orbit enumeration. Optimized.
  
- `MultiIndex` is an alternative implementation of `ind2sub` and `sub2ind`, designed to
  work with standard (and not cell) arrays. Optimized.

Organization
------------

License
-------
(C) Denis Rosset 2016
Placed under GPL license version 3, except for the files below, that can be distributed either
under the original license of the ported code, or the GPL version 3.

Code attribution:

- `LibDivide.m` is lifted from [libdivide](http://www.libdivide.com), under ZLIB

- `JacobiPolynomial.m` has parts taken from [QETLAB](http://www.qetlab.com/Main_Page), under 2-Clause BSD

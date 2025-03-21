# GSplines.jl

The package provide tools for the construction and use of Geometrically Smooth (G1) Splines .

## Example

Here is an example of construction of a G1 surface, and of a G1 basis associated to a mesh:

```
using GSplines, SemiAlgebraicTypes
m = offdata("cube.off")
gs = g1surface(m)
m1 = cc_subdivide(m)
bs = g1basis(m)
```

## Installation

It can be installed as follows:

```
] add https://github.com/AlgebraicGeometricModeling/GSplines.jl
```

[![](https://img.shields.io/badge/docs-blue.svg)](https://AlgebraicGeometricModeling.github.io/GSplines.jl/)
[![](https://img.shields.io/badge/source-orange)](https://github.com/AlgebraicGeometricModeling/GSplines.jl)

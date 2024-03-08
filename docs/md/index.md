# GSplines

The package provide tools for the construction and use of Geometrically Smooth (G1) Splines .


Here is an example of construction of a G1 surface, visualized with 'Axl':
```
using GSplines, SemiAlgebraicTypes, Axl

m = axldata("y1m1.axl")[1]  #read mesh from data file "y1m1.axl"
m[:color] = Axl.red
```

![y1m1](y1m1.png)

```
s = g1surface(m)
@axlview s
```

![y1g1](y1g1.png)

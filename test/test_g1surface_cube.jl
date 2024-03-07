using GSplines, SemiAlgebraicTypes, Axl

m  = offdata("cube.off")

#m  = offdata("1square.off")
#m  = offdata("2squares.off")

#hm = hmesh(m);
#cc_subdivide!(hm)
#divideEV(hm)

m[:color] = Axl.blue

s = g1surface(m)

@axlview s

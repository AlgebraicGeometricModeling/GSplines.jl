"""
 Script to test the construction of G1 basis functions
"""

using SemiAlgebraicTypes, GSplines

dir = joinpath(pwd(),"../data/off");

F = filter(x -> ( !startswith(x, "BEV_") & endswith(x, ".off") ) , readdir(dir))

#The input mesh must have EVs surrounded by just regular vertices; if this is not true, you must subdivide once (or more, for any reason).

for f in F[1:end]
    try
        local m = offread(joinpath(dir,f));
        m[:color] = Color(0,0,255);
        local hm = hmesh(m);
        divideEV(hm); #Subdivide mesh with Catmull-Clark subdivision scheme
        local t = @elapsed g1basis(m)
        @info "\033[96m$f\033[0m $t(s)"
    catch
        @warn "problem with $f"
    end
    println("\n")
end

#
# Script to test the G1 surface function
#

using SemiAlgebraicTypes, GSplines

dir = joinpath(G1S[:pkgdir],"data/axl");

F = filter(x -> endswith(x, ".axl") , readdir(dir))

# The input mesh must have EVs surrounded by just regular vertices; if this is not true, you must subdivide once (or more, for any reason).

S=["CS-S","CS-AS","NCS-S","NCS-AS"]; #Solving strategies

for f in F[1:end]
    try
        local m = axldata(f)[1];
        m[:color] = Color(0,0,255);
        local hm = hmesh(m);
        divideEV(hm); #Subdivide mesh with Catmull-Clark subdivision scheme
        for i in 1:4
            local t = @elapsed g1surface(hm, S[i]);
            local name =S[i]
            @info "\033[96m$f\033[0m $name    $t(s)"
        end
    catch
        @warn "problem with $f"
    end
    println("\n")
end

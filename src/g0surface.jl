export g0surface
"""
    g0surface(hm::HMesh,d=3), g0surface(m::Mesh,d=3)

This function takes in input a mesh (or an half egde data structure of a mesh) and returns a g0 multipatch surface obtained using the Approximate Catmull-Clark scheme (ACC3)

## Example

    using G1Splines
    m = offdata("cube.off")
    g0 = g0surface(m)

"""
function g0surface(hm::HMesh,d::Int=3)
    @assert d >= 3
    return acc_d(hm,d)
end

function g0surface(m::Mesh,d::Int=3)
    @assert d >= 3
    hm=hmesh(m);
    return acc_d(hm,d)
end

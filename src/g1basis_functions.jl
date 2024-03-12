export g1basis

"""
This function takes in input a quad mesh and returns a sparse matrix containing the coefficients defining a set of G1 biquintic basis functions on the input mesh.

 - ncols of the sparse matrix gives the dimension of the spline space,
 - nrows is the total number of control points in the mesh i.e. nfaces*36
"""

function g1basis(m::HMesh)
    basis=g1basis_bezier(m);
    return basis
end


function g1basis(m::Mesh)
    hm=hmesh(m)
    g1basis(hm);
end



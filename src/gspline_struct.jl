export GSpline
"""
Geometrically Smooth Splines

**Fields:**

    * degree: degree of the bspline functions on each patch
    * knots: vector of knots for each patch in each parameter direction
    * mesh: half-edge mesh defining the topological surface
    * ctrpoints: matrix of control points of size n*N, where
        - n is the number of coordinates of the control points (3 by default)
        - N is the number of control points (lsz * nbe where nbe is the number of half-edges of the mesh)
    * basis: array of sparse vectors representing the basis functions attached to the vertices and edges.
      The basis functions associated to the interior points of the faces are not represented.

"""

mutable struct GSpline
    regularity::Int64
    degree::Int64
    knots::Vector
    mesh::HMesh
    ctrpoints::Matrix{Float64}
#    basis::Vector{SparseVector{Float64,Int64}}
    glue::Dict{Int64,Float64}
    patchsize::Int64
    dimension::Int64
    interpolationglue::Dict{Int64,Array{Float64,1}}
    #interpolationglue::Dict{Int64,Float64}

    function GSpline()
        new(1, 0, Float64[], HMesh(),
            Matrix{Float64}(undef,0,0),
            Dict{Int64,Float64}(),
            0,0,
            Dict{Int64,Array{Float64,1}}()
            )
    end




"""
    Define the gspline structure where
        - r  is the regularity
        - kn is the sequence of knots on each face in the 2 directions
        - msh is the supporting defining the topology
        - P is the matrix of control points
    The degree d is computed from the number of times (d+1) the first knot is repeated.
"""
    function GSpline(r::Int64, kn::Vector, msh = HMesh(), P = Matrix{Float64}(undef,0,0))
        m, d = dim_deg(kn)

        gs = new(r, d, kn, msh,
                 P,
                 Dict{Int64,Float64}(),
                 m,
                 nbf(msh)*m*m,
                 Dict{Int64,Array{Float64,1}}()
                 )

        return gs
    end

end

export dim_deg
function dim_deg(kn)
    d = 1
    while kn[d] == kn[1]
        d+=1
    end
    d -= 2
    m  = length(kn)-d-1
    return m,d
end

export lsz, bsz, reg
# """
# Local size of an edge of GSpline.mesh.
# The number of functions attached to a hedge is lsz*lsz.
# """
function lsz(gs::GSpline)
    div(length(gs.knots)-gs.degree-1,2)
end

function bsz(gs::GSpline)
    length(gs.knots)-gs.degree-1
end


# """
# Regularity of the space of splines of gs.
# """
function reg(gs::GSpline)
    return gs.regularity
end


export bspline
"""
 Extract the bspline representation of face f
"""
function bspline(gs::GSpline, f::Int64, basis)
    
    m = bsz(gs)
    s = m*m
    points = fill(0.0, 3, m, m)
    e = face(gs.mesh,f)
    for i in 1:m
        for j in 1:m
            points[:, i, j] = gs.ctrpoints[:, idx(gs,e,i,j)]
        end
    end
    BSplineSurface(points, basis, basis)
end

"""
 Extract the bspline representations of the faces as a vector of bspline
 functions
"""
function bspline(gs::GSpline)
    bs = BSplineBasis(gs.knots,gs.degree+1,false)
    [bspline(gs,i, bs) for i in 1:nbf(gs.mesh)];
end



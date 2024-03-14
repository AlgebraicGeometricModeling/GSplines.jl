module GSplines

using LinearAlgebra
using SemiAlgebraicTypes
using Axl
using DataStructures
using ToeplitzMatrices
using FileIO
using SparseArrays

#currentdir =  @__DIR__

G1S = Dict{Symbol, Any}( :pkgdir => dirname(dirname(pathof(GSplines))) )

export G1S
    
#const MASKS=load(joinpath(currentdir,"EVsdict.jld2"))["EVsdict"];

const G1S[:masks] = Dict{String, Dict{Int64, Matrix{Float64}}}();

include("gspline_struct.jl")

include("indexing.jl")

include("g1constraint_matrix.jl")

include("g0surface.jl")
include("g1surface.jl")

include("acc_3.jl")
include("acc_d.jl")

include("ordering.jl")

include("mesh_operations.jl")
include("mesh_boundary.jl")

include("mask_degree_elevate.jl")
include("mask_operations.jl")
include("order_bezier_cp.jl")

include("g1basis.jl")
include("g1basis_bezier.jl")

include("get_basis_from_sparse.jl")

include("g1dimension.jl")

include("eval.jl")
include("g1basis_assemble_spline.jl")

include("insert_knot.jl")

include("g1basis_gluingdata.jl")
include("g1basis_assemble_gluing.jl")



include("off_read.jl")
include("gismo_write.jl")
include("axl_write.jl")

end

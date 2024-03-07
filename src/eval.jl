export eval_vec, eval_bs

# """
# Compute the vector of values of bspline functions at the point of face f with parameters u=[u1, u2]
# They are numbered according to their index idx(gs, e, i, j) where e is the
# first half-edge of the face.
# """
function eval_vec(gs::GSpline, f::Int64, u::Vector )
    bs = BSplineBasis(gs.knots,gs.degree+1,false)
    b1, r1 = eval_rng(bs,u[1])
    b2, r2 = eval_rng(bs,u[2])
    e = gs.mesh.faces[f]

    vct = spzeros(size(gs.ctrpoints,2))
    for (i,s) in zip(r1,b1), (j,t) in zip(r2, b2)
        vct[idx(gs, e, i, j)] = s*t
    end
    vct
end


function eval_mat(bs::BSplineBasis, u::Vector, d1=0, d2=0 )

    b1, r1 = eval_rng(bs,u[1],d1)
    b2, r2 = eval_rng(bs,u[2],d2)

    M = zeros(typeof(b1[1]),length(bs), length(bs))
    #vct = spzeros(size(gs.ctrpoints,2))
    for (i,s) in zip(r1,b1), (j,t) in zip(r2, b2)
        M[i,j]= s*t
    end
    return M
end



# """
# Compute the value of gs at the point of face f with parameters u=[u1, u2]
# """
function Core.eval(gs::GSpline, f::Int64, u::Vector)
    gs.ctrpoints*eval_vec(gs, f, u)
end


# """
# Compute the value of the basis elements Bs at the point of face f with parameters u=[u1, u2]
# """
function eval_bs(Bs, gs::GSpline, f::Int64, u::Vector)
    Bs'*eval_vec(gs,f,u)
end



# """
# Compute the vector of values of all the g1 basis elements at point (f,u)
#   - the first element of the basis are gs.basis
#   - the last elements are the inner face bsplines, numbered by blocks of
#     (2l-4)^2, one for each face.
# """
function eval_basis(gs::GSpline, f, u)

    reg = gs.regularity
    v = eval_vec(gs, f, u)
    r = spzeros(dim(gs))
    for i in 1:length(gs.basis)
        r[i] = dot(v,gs.basis[i])
    end

    l = lsz(gs)
    s = length(gs.basis) + 1 + (f-1)*4*(l-reg-1)^2

    e = gs.mesh.faces[f]
    for j in reg+2:2*l-reg-1
        for k in reg+2:2*l-reg-1
            r[s] = v[idx(gs, e, j, k)]
            s = s+1
        end
    end
    r
end

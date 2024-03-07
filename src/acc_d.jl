export acc3,acc_d

function face_point(msh, e, v)
    p =  point_of(msh,e)*v
    p += point_of(msh,next(msh,e))*2
    p += point_of(msh,prev(msh,e))*2
    p += point_of(msh,next(msh,next(msh,e)))
    p /= (v+5)
    return p
end


function acc3(gs::GSpline)
    msh = gs.mesh
    l = lsz(gs)
    s = l*l
    N = nbe(msh)*s
    C  = fill(0.0, 3, N)

    v = 0

    CCE =  SemiAlgebraicTypes.ccw_edges(msh)
    msh.ccw_e = CCE

    for p in 1:nbv(msh)
        if length(CCE[p])>0

            E = CCE[p]
            v = length(E)
            if opp(msh,E[1]) != 0
                # Interior vertex
                p0 = msh.points[:,p]*v*v
                for e in E
                    p0 += point_of(msh, next(msh,e))*4
                    p0 += point_of(msh, next(msh, next(msh,e)))
                end
                p0 /= v*(v+5)
            elseif opp(msh,prev(msh,E[1])) != 0
                # Boundary non-corner vertex
                p0 = msh.points[:,p]*4
                p0 += point_of(msh, next(msh,E[1]))
                p0 += point_of(msh, prev(msh,E[end]))
                p0 /= 6
            else
                # Corner vertex
                p0 = msh.points[:,p]
            end
            for e in E
                C[:,idx(gs, e,1,1)] = p0
            end
        end
    end

    for e in 1:nbe(msh)
        # Edge point
        if opp(msh,e) != 0
            if opp(gs.mesh,CCE[ptidx_of(msh,e)][1])!=0
                v = length(CCE[ptidx_of(msh,e)])
            else
                v = length(CCE[ptidx_of(msh,e)])*2
            end
            p1  = point_of(msh,e)*2*v
            p1 += point_of(msh, next(msh,e))*4
            p1 += point_of(msh, prev(msh,e))*2
            p1 += point_of(msh,next(msh,next(msh,e)))

            no = next(msh,opp(msh,e))
            p1 += point_of(msh,next(msh,no))*2
            p1 += point_of(msh,next(msh,next(msh,no)))
            p1 /= (v+5)*2
            C[:,idx(gs, no,1,2)] = p1
            C[:,idx(gs, e, 2,1)] = p1

        else

            v = 4
            p1  = point_of(msh,e)*2
            p1 += point_of(msh, next(msh,e))
            p1 /= 3
            C[:,idx(gs, e,2,1)] = p1
            p1  = point_of(msh, next(msh,e))*2
            p1 += point_of(msh,e)
            p1 /= 3
            C[:,idx(gs, e,3,1)] = p1
        end

        # Face point
        C[:,idx(gs, e, 2,2)] = face_point(msh,e,v)

    end
    return C
end

function acc3(hm::HMesh)
    knts = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    gs = GSpline(0, knts, hm )
    gs.ctrpoints = acc3(gs)
    return gs
end




function acc_d(msh::HMesh,d::Int64)
    if d<3
        return
    end
    gs0 = acc3(msh)
    m=d+1
    knts = zeros(1,2*(m))
    knts[1:m].=0
    knts[m+1:end].=1

    gs = GSpline(0,vec(knts),msh)

    N = nbf(msh)*m*m
    gs.ctrpoints = fill(0.0, 3, N)

    M = fill(zero(gs0.ctrpoints[:,1]),4,4)
    for f in 1:nbf(msh)
        ef = face(msh,f)
        for i in 1:4
            for j in 1:4
                M[i,j] = gs0[ef,i,j]
            end
        end

        M0=copy(M)
        M1=[];
        push!(M1,M0)
        for i in 1:d-3
            push!(M1,degree_elevate(copy(M1[end])))
        end

        for i in 1:m
            for j in 1:m
                gs[ef,i,j] = M1[end][i,j]
            end
        end


    end
    return gs
end

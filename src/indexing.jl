export idx, oidx, is_g1_interior
# """
# Index of the bspline function associated to the edge e and position (i, j)
# with 1 <= i,j <= 2*lsz
# """
function oidx(gs, e, i, j)
    l = lsz(gs)
    if i<= l
        if j<= l
            return (e-1)*l*l + (j-1)*l+i
        else
            p = prev(gs.mesh,e)
            return (p-1)*l*l + (i-1)*l+(2*l-j+1)
        end
    else
        if j<= l
            n = next(gs.mesh,e)
            return (n-1)*l*l + (2*l-i)*l+j            
        else
            n = next(gs.mesh,next(gs.mesh,e))
            return (n-1)*l*l + (2*l-j)*l+(2*l-i+1)
        end
    end
end

function idx(gs::GSpline, e, i, j)
    m  = bsz(gs)
    f  = hedge(gs.mesh, e).face
    ef = face(gs.mesh, f)

    if ef == e
        
        return (f-1)*m*m + (j-1)*m + i
        
    elseif prev(gs.mesh,e) == ef
        
        return (f-1)*m*m + (i-1)*m + (m-j+1)

    elseif next(gs.mesh,e) == ef
        
        return (f-1)*m*m + (m-i)*m + j            

    else
        
        return (f-1)*m*m + (m-j)*m + (m-i+1)

    end
end

function idx(msh, m, e, i, j)
    f  = hedge(msh, e).face
    ef = face(msh, f)

    if ef == e
        
        return (f-1)*m*m + (j-1)*m + i
        
    elseif prev(msh,e) == ef
        
        return (f-1)*m*m + (i-1)*m + (m-j+1)

    elseif next(msh,e) == ef
        
        return (f-1)*m*m + (m-i)*m + j            

    else
        
        return (f-1)*m*m + (m-j)*m + (m-i+1)

    end
end


function fidx(gs, corners, n)
    m  = bsz(gs)

    f = div(n-1,m*m)+1
    corner = corners[f]
    
    if corner == 1
        
        i = rem(n-1,m)+1
        j = div(rem(n-1,m*m),m)+1
        
    elseif corner == 2

        j = m-rem(n-1,m)
        i = div(rem(n-1,m*m),m)+1

    elseif  corner == 3

        j = rem(n-1,m)+1
        i = m-div(rem(n-1,m*m),m)

    else
        i = m-rem(n-1,m)
        j = m-div(rem(n-1,m*m),m)
    end

    return f,i,j
end

function lidx(b, e, i, j=1)
    return (e-1)*b*b + (j-1)*b+i
end

function bidx(l, i, j)
    if i==1
        return [(j-1)%l+1,div((j-1),l)+1]
    elseif i==2
        return [2*l-div((j-1),l),(j-1)%l+1]
    elseif i==3
        return [2*l-(j-1)%l,2*l-div((j-1),l)]
    else
        return [div((j-1),l)+1,2*l-(j-1)%l]
    end
end


#=
export idx

function idx(gs::GSpline, e, i, j)
    m  = bsz(gs)
    f  = hedge(gs.mesh, e).face
    ef = face(gs.mesh, f)

    if ef == e
        
        return (f-1)*m*m + (j-1)*m + i
        
    elseif prev(gs.mesh,e) == ef
        
        return (f-1)*m*m + (i-1)*m + (m-j+1)

    elseif next(gs.mesh,e) == ef
        
        return (f-1)*m*m + (m-i)*m + j            

    else
        
        return (f-1)*m*m + (m-j)*m + (m-i+1)

    end
end

function idx(msh, m, e, i, j)
    f  = hedge(msh, e).face
    ef = face(msh, f)

    if ef == e
        
        return (f-1)*m*m + (j-1)*m + i
        
    elseif prev(msh,e) == ef
        
        return (f-1)*m*m + (i-1)*m + (m-j+1)

    elseif next(msh,e) == ef
        
        return (f-1)*m*m + (m-i)*m + j            

    else
        
        return (f-1)*m*m + (m-j)*m + (m-i+1)

    end
end
=#


# """
# Return true iff the index i of the bspline basis function is 
# interior to a face.
# """
function is_g1_interior(gs, i)
    l = lsz(gs)
    s = l*l
    n = (i-1)%s+1
    a = (n-1)%l+1
    b = (n-1)÷l+1
    return a>2 && b>2 
end

function set_interior!(gs, C1, C0, d=1)
    m = bsz(gs)
    for f in 1:nbf(gs.mesh)
        for i in d+2:m-d-1
            for j in d+2:m-d-1
                C1[:, (f-1)*m*m+(j-1)*m+i] = C0[:,(f-1)*m*m+(j-1)*m+i]
            end
        end
    end
end

function Base.getindex(gs::GSpline,e,i,j)
 gs.ctrpoints[:,idx(gs,e,i,j)]
end

function Base.setindex!(gs::GSpline, v, e,i,j)
    gs.ctrpoints[:,idx(gs,e,i,j)] = v
end

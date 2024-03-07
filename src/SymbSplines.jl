
function wof(t, i, j, var)
    if t[i+1] != t[j+i+1]
        return (var-t[i+1])*(1/(t[i+1+j]-t[i+1]))
    else
        return 0
    end
end

function search_in(v, tau)
    local i;
    for i in 1:length(tau)
        if v == tau[i]
            return i
        end
    end
    println("Knot ", v," don't exists in ",tau)
end
    
function associated_subdivision(t)
    local tau, cmp, i;
    tau = [];
    cmp = 1;
    tau = [tau; t[1]];
    for i in 1:length(t)
        if tau[cmp] != t[i]
            cmp = cmp+1;
            tau = [tau; t[i]]
        end
    end
    return tau
end


function bspline_rec(i, d, kn, var)
    local tau, B0, index;
    if length(kn)-d-1 < i || i < 0
        println("incompatible indexes");
        return;
    end
    if d == 0 
        tau = associated_subdivision(kn);
        B0 = fill(0, length(tau)-1);
        if kn[i+1] != kn[i+2]
            index = search_in(kn[i+1], tau);
            B0[index] = 1
        end
        return B0
    end
    return wof(kn, i, d, var)*bspline_rec(i, d-1, kn, var)+(1-wof(kn, i+1, d, var))*bspline_rec(i+1, d-1, kn, var)
end

function bspline(i, kn, u)
    m, d = dim_deg(kn)
    bspline_rec(i,d,kn,u)
end

function bspline(i, j, kx, ky, u, v)

    m1,d1 = dim_deg(kx)
    m2,d2 = dim_deg(ky)
    
    Bx = bspline(i, kx, u);
    By = bspline(j, ky, v);
    BB = fill(zero(Bx[1]), length(Bx), length(By));
    for l in 1:length(Bx) 
        for m in 1:length(By) 
            BB[l, m] = Bx[l]*By[m]
        end
    end
    return BB
end

function generic_bspline(k1, k2, u, v, A)

    m1, d1  = dim_deg(k1)
    m2, d2  = dim_deg(k2)
    sum(A[i,j]*bspline(i-1,j-1,k1,k2,u,v)  for i in 1:m1, j in 1:m2)
    
end
    
   
function a_glue(a0, a1, u)
    a0*(u-1)^2 - a1*u^2
end

function coefficients(E, u)

    e = typeof(E)[]
    e=fill(zero(E), maxdegree(E,u)+1)
    for t in E
        e[maxdegree(t,u)+1] += subs(t,u=>1)
    end
    e
end

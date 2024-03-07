#include("../src/G1Splines.jl")
#include("../src/eval.jl")


function Aglue(N1,N2,x)
    a = 2*cos(2*pi/N1)*(1-x)^2-2*cos(2*pi/N2)*x^2 #+ x^2
end

function AAgluex(N1,N2,x)
    2*x-1 #a = 2*cos(2*pi/N1)*(1-x)^2-2*cos(2*pi/N2)*x^2 #+ x^2
end


function Agluex(N1,N2,x)
    if N1==3
        c1 =-1//2
    elseif N1 == 4
        c1 = 0//1
    else
        #@warn "glueing data non-exact at start point"
        c1=cos(2*pi/N1)
    end
    if N2==3
        c2 =-1//2
    elseif N2 == 4
        c2 = 0//1
    else
        #@warn "glueing data non-exact at end point"
        c2=cos(2*pi/N2)
    end
    a = 2*c1*(1-x)^2-2*c2*x^2 #+ x^2
#    a = 2*c1*(1-x)-2*c2*x #+ x^2
end


function f2_idx(m,k,i,j)
    return (k-1)*m^2+(j-1)*m+i
end


"""
Compute the matrix defining the edge spline basis, in a compressed form.
"""
function g1matrix_edge(kn::Array, N1::Int64, N2::Int64;
                       zeroOnEdge=false,
                       equalAcrossEdge=false,
                       oppositeAcrossEdge=false,
                       C1VirtualEdge=false)

    m, deg = dim_deg(kn)
    println("valences: ", N1, " ", N2)


    println("m = ", m, " d = ", deg)

    T = typeof(kn[1])

    bs = BSplineBasis(kn,deg+1,false)


    #U = [sum(kn[i:i+deg])//(deg) for i in 1:length(kn)-deg]
    #s = div(5*length(kn),4)+5
    s = length(kn)+5
    U = [i//s for i in 0:s]

    #println(U)

    N = 2
    P = N*m*m

    L = length(U) + m + 2*m*(m-2) + 18 + (m-4)


    I = Int64[]

    S = spzeros(T,L,P)

    l = 1

    # Equality of coefficients on the shared edges ---> C0 conditions
    for i in 1:m
            S[l,f2_idx(m,2,i,1)] =  1
            S[l,f2_idx(m,1,1,i)] = -1
            l+=1
            #println(": ",idx(gs,e,i,1),"=",idx(gs,n,1,i))
    end

    b = -1

    for u in U

        a = Agluex(N1,N2,u)
        #println(a)

        M1 = eval_mat(bs, [u, zero(T)], 0,1)
        M2 = b*eval_mat(bs, [zero(T),u], 1,0) + a*eval_mat(bs, [zero(T),u], 0, 1)

        #println(typeof(M1), "\n", typeof(M2),"\n\n")

        for i in 1:m
            for j in 1:m
                if M1[i,j] != 0
                    #S[l, idx(hm,m,e,i,j)] = - M1[i,j]
                    S[l, f2_idx(m,2,i,j)] = - M1[i,j]
                end
                if M2[j,i] != 0
                    #S[l, idx(hm,m,n,j,i)] = M2[j,i]
                    S[l, f2_idx(m,1,j,i)] = M2[j,i]
                end
            end
        end
        l+=1
    end


    #Zero values for the inner face ctpoints
    if true
        for i in 1:m
            for j in 3:m
                push!(I, f2_idx(m,2,i,j))
                push!(I, f2_idx(m,1,j,i))
            end
        end
    end


    # Equal pair of ctpoints across the edge for higher derivatives (before the crossing vertex relations)
    if equalAcrossEdge
        for i in 3:m-2
            S[l, f2_idx(m,2,i,2)]= 1
            S[l, f2_idx(m,1,2,i)]= -1
            l+=1
            #=if i%5==1 #C1 contraints across the virtual edge
                p=i;
                S[l, f2_idx(m,1,2,p-1)]= 1
                S[l, f2_idx(m,1,2,p+1)]= 1
                S[l, f2_idx(m,1,2,p)]= -2
                l+=1
                S[l, f2_idx(m,1,1,p-1)]= 1
                S[l, f2_idx(m,1,1,p+1)]= 1
                S[l, f2_idx(m,1,1,p)]= -2
                l+=1
                S[l, f2_idx(m,2,p-1,2)]= 1
                S[l, f2_idx(m,2,p+1,2)]= 1
                S[l, f2_idx(m,2,p,2)]= -2
                l+=1

                #b30=0 for first ker basis
                    #S[l, f2_idx(m,2,p-3,1)]= 1
                    #l+=1
                #S[l, f2_idx(m,1,1,p+3)]= 1
                #l+=1
                #S[l, f2_idx(m,2,p,1)]= 1
                #l+=1
                    #S[l, f2_idx(m,1,1,p)]= 1
                    #l+=1
                #=S[l, f2_idx(m,2,p-3,1)]= 1
                #l+=1
                S[l, f2_idx(m,2,p+3,1)]= 1
                l+=1=#


            end=#
        end
    end



    # Oposite pair of ctpoints across the edge for higher derivatives (before the crossing vertex relations)
    if oppositeAcrossEdge
        for i in 3:m-2
            S[l, f2_idx(m,2,i,2)]= 1
            S[l, f2_idx(m,1,2,i)]= 1
            l+=1
        end
    end


    if C1VirtualEdge
        for i in 3:m-2

            if i%5==1 #C1 contraints across the virtual edge
                p=i;
                S[l, f2_idx(m,1,2,p-1)]= 1
                S[l, f2_idx(m,1,2,p+1)]= 1
                S[l, f2_idx(m,1,2,p)]= -2
                l+=1
                S[l, f2_idx(m,1,1,p-1)]= 1
                S[l, f2_idx(m,1,1,p+1)]= 1
                S[l, f2_idx(m,1,1,p)]= -2
                l+=1
                S[l, f2_idx(m,2,p-1,2)]= 1
                S[l, f2_idx(m,2,p+1,2)]= 1
                S[l, f2_idx(m,2,p,2)]= -2
                l+=1
            end
        end
    end

    # Zero for the edge coefficients
    if zeroOnEdge
        for i in 1:m
            S[l, f2_idx(m,2,i,1)]= 1
            l+=1
            S[l, f2_idx(m,1,1,i)]= 1
            l+=1
        end
    end

    # 8 zero ctps at end vertex
    if true
        for i in 1:2
            for j in m-1:m
                push!(I, f2_idx(m,1,i,j))
            end
        end
        for i in 1:2
            for j in m-1:m
                push!(I, f2_idx(m,2,j,i))
            end
        end
    end


    # 8 zero ctps at first vertex
    if true
        for k in 1:2
            for i in 1:2
                for j in 1:2
                    push!(I, f2_idx(m,k,i,j))
                end
            end
        end
    end


    #println(sort(I))

    J = deleteat!([i for i in 1:P], sort(I))

    println("Size ", size(S), "   reduced: ", l, "  ", length(J))

    #println("DIMENSION: " , P - rankx(Matrix(S)))

    S[1:l,J], J
end

function expand_spmatrix(K,I,m)
    Ke = spzeros(2*m^2,size(K,2))
    for j in 1:size(K,2)
        for i in 1:size(K,1)
            if K[i,j] != 0
                Ke[I[i],j] = K[i,j]
            end
        end
    end

    return Ke
end

"""
 Compute the matrix defining the G1 splines on the mesh hm with the knot distribution kn using quadratic glueing data.

 It assumes that boundary vertices are regular.
"""
function g1matrix(hm, kn)
    m, d = dim_deg(kn)

    bs = BSplineBasis(kn,d+1,false)

    #gs = GSpline(0, kn, hm )
    T = typeof(kn[1])

    s = m+7
    U = [i//s for i in 0:s]

    N = nbf(hm)*m*m

    ne = 0
    for e in 1:nbe(hm)
        o = opp(hm,e)
        if e < o  #pair of shared edge with lower first index
            ne+=1
        end
    end

    L =ne*(length(U)+m)
    #L =ne*(m)

    E = ccw_edges(hm)

    S = spzeros(T,L,N)
    l = 1
    # Equality of coefficients on the shared edges
    for e in 1:nbe(hm)
        o = opp(hm,e)
        if e < o  #pair of shared edge with lower first index
            n = next(hm,o)
            for i in 1:m
                S[l,idx(hm,m,e,i,1)] = 1
                S[l,idx(hm,m,n,1,i)] = -1
                l+=1
                #println(": ",idx(gs,e,i,1),"=",idx(gs,n,1,i))
            end
        end
    end
    #println("rank S:", rank(S))

    b = -1
    for e in 1:nbe(hm)
        o = opp(hm,e)
        if e < o  #pair of shared edge with lower first index
            N1 = length(E[edge(hm,e).point])
            if N1 == 2  N1=4 end
            N2 = length(E[edge(hm,o).point])
            if N2 == 2  N2=4 end

            #A = Aglue(V1,V2)
            #        A = -(x-1//1)^2

            #println("edge valences: ",N1, " ",N2)

            n = next(hm,o)


            for u in U

                a = Agluex(N1,N2,u)

                M1 = eval_mat(bs, [u, 0//1], 0,1)
                M2 = b*eval_mat(bs, [0//1,u], 1,0) + a*eval_mat(bs, [0//1,u], 0,1)

                #println(l,"\n",M1,"\n")
                #println(M2)
                #println("\n--------------------------------------")
                #println(M2)
                for i in 1:m
                    for j in 1:m
                        if M1[i,j] != 0
                            S[l, idx(hm,m,e,i,j)] = - M1[i,j]
                        end
                        if M2[j,i] != 0
                            S[l, idx(hm,m,n,j,i)] = M2[j,i]
                        end
                    end
                end
                l+=1
            end
        end
    end

    #r = rankx(Matrix(S));
    #println(size(S), " ",svd(Matrix(S)).S)
    #dm = size(S,2)-r;

    return S;
end

function g1dimension(hm,kn)
    S = g1matrix(hm,kn)
    return size(S,2)-rank(S)
end

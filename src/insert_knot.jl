#=
function knots(d,N,mu=1)
    kn = [0.0]

    for i in 1:d push!(kn,0) end
    for i in 1:N-1
        for k in 1:mu
            push!(kn,i/N)
        end
    end
    for i in 1:(d+1) push!(kn,1) end
    return kn
end
=#

function insert_knot(M::Array, knt::Array, t::Float64)

    m, d = dim_deg(knt);
    nbelem=size(M,1);
    nbrows=size(M,2);
    nbcols=size(M,3);
    Mu=zeros(nbelem,nbrows,nbcols+1);
    Mv=zeros(nbelem,nbrows+1,nbcols+1);

    k=1
    while !(knt[k]<=t<knt[k+1])
        k+=1
    end

    knt_new=zeros(length(knt)+1);
    knt_new[1:k]=knt[1:k];
    knt_new[k+1]=t;
    knt_new[k+2:end]=knt[k+1:end];

    #Compute the alphas
    alpha=[];

    for i in k-d+1:k
        push!(alpha,((t-knt[i])/(knt[i+d]-knt[i])));
    end

    #Knot insertion in u dir
    for j in 1:nbrows
        l=1;
        for h in 1:nbcols+1
            if h<=k-d
                Mu[:,j,h]=M[:,j,h];
            elseif h>=k+1
                Mu[:,j,h]=M[:,j,h-1];
            else
                Mu[:,j,h]=alpha[l]*M[:,j,h]+(1-alpha[l])*M[:,j,h-1];
                l+=1;
            end
        end
    end


    #Knot insertion in v dir
    for j in 1:nbrows+1
        l=1;
        for h in 1:nbcols+1
            if h<=k-d
                Mv[:,h,j]=Mu[:,h,j];
            elseif h>=k+1
                Mv[:,h,j]=Mu[:,h-1,j];
            else
                Mv[:,h,j]=alpha[l]*Mu[:,h,j]+(1-alpha[l])*Mu[:,h-1,j];
                l+=1;
            end
        end
    end



    return Mv,knt_new



end


#=function deCasteljau(ctrpts::Array,u0::Float64)

#ctrpts is a 3xn matrix

n=size(ctrpts)[2];
M=zeros(n,n,3);
M[:,1,1]=ctrpts[1,:];
M[:,1,2]=ctrpts[2,:];
M[:,1,3]=ctrpts[3,:];

for i in 2:n
    for j in 2:n
        M[j,i,:]=(1-u0)*M[j-1,i-1,:]+u0*M[j,i-1,:];
    end
end

new1=diag(M[:,:,:]);
new2=reverse(M[end,:,:]);

return new1,new2

end=#



function deCasteljau(ctrpts::Array,u0::Float64)

#ctrpts is a 1xn matrix

n=length(ctrpts);
M=zeros(n,n);
M[:,1]=ctrpts;

for i in 2:n
    for j in 2:n
        M[j,i]=(1-u0)*M[j-1,i-1]+u0*M[j,i-1];
    end
end

new1=diag(M);
new2=reverse(M[end,:]);

return new1,new2

end

function glob_gluing_subd(ctrpts::Array, nbsubd::Int64)

M=fill(Vector{Any}(), nbsubd+1, 2^nbsubd);
u0=0.5;

M[1,1]=ctrpts;

if nbsubd==0
    return M
else

    for i in 1:nbsubd
        for j in 1:2^i
            idxglu=(j+1)%2+1;
            M[i+1,j]=deCasteljau(M[i,div(j+1,2)],u0)[idxglu];
        end
    end

end

return M[end,:]

end


function transform_knot_vect(M::Array, knt_in::Array, knt_fin::Array)

    if knt_fin==knt_in
        return M
    else

        m,d=dim_deg(knt_in);

        knts=knt_fin[(d+2):(end-(d+1))];

        for i in 1:length(knts)

            sol=insert_knot(M,knt_in,knts[i]);

            M=sol[1];
            knt_in=sol[2];

        end

    end


    return M;

end


function bezier_to_spline(B, knt::Array, N::Int64, S::String)


    nbbas=size(B,2);
    m,d=dim_deg(knt);
    knt_bez=knots(d,1,1);
    m_bez,d=dim_deg(knt_bez);
    newbasis=zeros(N*m*m,nbbas);
    loc_bas=zeros(1,m_bez,m_bez);

    for k in 1:nbbas
        bas=B[:,k];
        tmp_mat=zeros(N,m_bez,m_bez);
        data=findnz(bas);
        row=data[1];
        val=data[2];

        for i in 1:length(row)
            j=div(row[i],m_bez^2)+1
            h=div(row[i]%m_bez^2,m_bez)+1;
            w=row[i]%m_bez;
            tmp_mat[j,w,h]=val[i];
        end


        for l in 1:N
            loc_bas[1,:,:]=tmp_mat[l,:,:];
            nCP=transform_knot_vect(loc_bas, knt_bez, knt);


            for r in 1:m
                for c in 1:m
                    if false
                        newbasis[loc_idx(m,l,r,c),k]=nCP[1,r,c];
                    else
                        if S=="innervertex"
                            if c<3 || r<3
                                newbasis[loc_idx(m,l,r,c),k]=nCP[1,r,c];
                            else
                                newbasis[loc_idx(m,l,r,c),k]=0;
                            end
                        elseif S=="bordervertex"
                            if l==1
                                if r<3
                                    newbasis[loc_idx(m,l,r,c),k]=nCP[1,r,c];
                                else
                                    newbasis[loc_idx(m,l,r,c),k]=0;
                                end
                            elseif l==2
                                if c<3
                                    newbasis[loc_idx(m,l,r,c),k]=nCP[1,r,c];
                                else
                                    newbasis[loc_idx(m,l,r,c),k]=0;
                                end
                            end
                        end
                    end
                end
            end


        end

    end

    return sparse(newbasis);

end

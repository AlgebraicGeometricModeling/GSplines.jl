export degree_elevate, degree_elevate_mask 

degree_elevate =function(M)
    MM = fill(zero(M[1,1]), size(M,1)+1, size(M,2)+1)
    n = size(M,2)-1
    m = size(M,2)-1
    for i in 1:n+2
        for j in 1:m+2
            if ((i-1)*(j-1))!=0
                MM[i,j]+=((i-1)*(j-1)/((n+1)*(m+1)))*M[i-1,j-1]
            end
            if ((n+1-(i-1))*(j-1))!=0
                MM[i,j]+=((1-(i-1)/(n+1))*((j-1)/(m+1)))*M[i,j-1]
            end
            if (i-1)*(m+1-(j-1))!=0
                MM[i,j]+=(((i-1)/(n+1)))*(1-(j-1)/(m+1))*M[i-1,j]
            end
            if ((n+1-(i-1))*(m+1-(j-1)))!=0
                MM[i,j]+=((1-(i-1)/(n+1))*(1-(j-1)/(m+1)))*M[i,j]
            end
        end
    end
    return MM
end

function degree_elevate_mask(N::Int, t::String, k::Int=3)
    #t::type of the patch, t==1 interior t==2 border t==3 rightEVpatch t==4 leftEVpatch t==5 innerface with borderEV t==6 corner
    #k::border vertex valence
    B=acc3bigmatrix(N,"INNEREV");
    
    if t=="INNEREV" #1
        B=acc3bigmatrix(N,"INNEREV");
    elseif t=="BORDER1" #2
        B=acc3bigmatrix(N,"BORDER1");
    elseif t=="BORDEREVR" #3
        B=acc3bigmatrix(N,"BORDEREVR");
    elseif t=="BORDEREVL" #4
        B=acc3bigmatrix(N,"BORDEREVL");
    elseif t=="BORDERFACEEV" #5
        B=acc3bigmatrix(N,"BORDERFACEEV");
    elseif t=="CORNER" #6
        B=acc3bigmatrix(N,"CORNER");
    elseif t=="BORDER2" #7
        B=acc3bigmatrix(N,"BORDER2");
    else
        @warn "Type not recognized, using default one"
    end

    D=fill(Vector{Float64}(),4,4);
    D[1,1]=B[:,1];
    D[1,2]=B[:,2];
    D[1,3]=B[:,7];
    D[1,4]=B[:,5];
    D[2,1]=B[:,3];
    D[2,2]=B[:,4];
    D[2,3]=B[:,8];
    D[2,4]=B[:,6];
    D[3,1]=B[:,14];
    D[3,2]=B[:,16];
    D[3,3]=B[:,12];
    D[3,4]=B[:,11];
    D[4,1]=B[:,13];
    D[4,2]=B[:,15];
    D[4,3]=B[:,10];
    D[4,4]=B[:,9];

    Q1=degree_elevate(D);
    Q2=degree_elevate(Q1); #biquintic masks points

return Q2;

end

"""
  Elevate the mask of degree 3x3 to degree 5x5
"""
function degree_elevate_mask(B::Matrix)

    D=fill(Vector{Float64}(),4,4);
    D[1,1]=B[:,1];
    D[1,2]=B[:,2];
    D[1,3]=B[:,7];
    D[1,4]=B[:,5];
    D[2,1]=B[:,3];
    D[2,2]=B[:,4];
    D[2,3]=B[:,8];
    D[2,4]=B[:,6];
    D[3,1]=B[:,14];
    D[3,2]=B[:,16];
    D[3,3]=B[:,12];
    D[3,4]=B[:,11];
    D[4,1]=B[:,13];
    D[4,2]=B[:,15];
    D[4,3]=B[:,10];
    D[4,4]=B[:,9];

    Q1=degree_elevate(D);
    Q2=degree_elevate(Q1); #biquintic masks points

    return Q2;

end

export g1surface
"""
    g1surface(hm::HMesh,S::String="CS-S")

This function takes in input a mesh in the half egde data structure and a string with the
G1 solving strategy and will return an array composed by the patches constituing the surface.

By default, if two EVs are connected, the mesh is split.
If the option `check_ev = false`, this does not happen, but the construction may be wrong.

The input string contains the solving strategy for the construction of the G1 surface to be selected from the
    following four: "CS-S","CS-AS","NCS-S","NCS-AS". The default strategy is "CS-S".

## Example

    using G1Splines
    m = offdata("cube.off")
    g1 = g1surface(m)
    
"""
function g1surface(hm::HMesh, S::String = "CS-S"; check_ev = true)

    if nbf(hm)<3
        return acc_d(hm,5)
    end

    if check_ev
        divideEV(hm)
    end
    
    knts5 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    gs5 = GSpline(0, knts5, hm )

    ncp_tot=36*nbf(hm);
    C5  = fill(0.0, 3, ncp_tot) #global matrix of cp

    P = neigh(hm);

    println(" P = ", P)
    
    CP3=zeros(3,16,1);
    CP5=zeros(3,36,1);
    #sup=[]; #not needed ?

    D=collect(keys(get(G1S[:masks],S,0)));

    #Inner patches with, eventually, EVs not on boundaries

    for i in 1:length(P[1])
        N=P[1][i][1]; #EV valence
        edge=P[1][i][2]; #face hedge going out from the EV (when present)
        if N==4
            #if D[1]==1
            #    B=acc3bigmatrix(N,"INNEREV")
            #else
            B=get(get(G1S[:masks],"RV",Dict{Int,Any}()),N,acc3bigmatrix(N,"INNEREV")) #Bézier ACC3 subdivision matrix
            mpoints=hm.points[:,P[1][i][3:end]] #Mesh points corresponding to the neighborhood
            CP3[:,:,1]=mpoints*B; #Bézier control points
            order_bezier_cp(gs5,C5,mpoints*B,edge,3);
            #push!(sup,order_plot(CP3,3));
            #end
        else
            #elseif N in D
            B=get(get(G1S[:masks],S,Dict{Int,Any}()),N,G1matrix(N,S,"INNER"));
            #else
            #    B=G1matrix(N,S,"INNER");
            #end
            mpoints=hm.points[:,P[1][i][3:end]]; #Mesh points corresponding to the neighborhood
            CP5[:,:,1]=mpoints*B; #bezier control points
            order_bezier_cp(gs5,C5,mpoints*B,edge,5);
            #push!(sup,order_plot(CP5,5));
        end
    end

    #Regular border patches

    for i in 1:length(P[2])
        N=P[2][i][1]; #First vertex valence regular vertex
        edge=P[2][i][2]; #face hedge going out from the EV (when present)
        #k=2; #Two faces on each adjacent vertex
        #if N==4
        B=get(get(G1S[:masks],"BP",Dict{Int,Any}()),N,acc3bigmatrix(N,"BORDER1"));
        mpoints=hm.points[:,P[2][i][3:end]];
        CP3[:,:,1]=mpoints*B;
        order_bezier_cp(gs5,C5,mpoints*B,edge,3);
        #order_bezier_cp(gs3,C3,mpoints*B,edge,3);
        #push!(sup,order_plot(CP3,3));
    #=elseif N!=4 && opp(hm,next(hm,next(hm,edge)))==0
        B=G1matrix(N,S,"EVREGBORDER1");
        mpoints=hm.points[:,P[2][i][3:end]]; #Mesh points corresponding to the neighborhood
        order_bezier_cp(gs5,C5,mpoints*B,edge,5);
        CP5[:,:,1]=mpoints*B; #bezier control points
        push!(sup,order_plot(CP5,5));
    elseif N!=4 && opp(hm,next(hm,edge))==0
        B=G1matrix(N,S,"EVREGBORDER2");
        mpoints=hm.points[:,P[2][i][3:end]]; #Mesh points corresponding to the neighborhood
        order_bezier_cp(gs5,C5,mpoints*B,edge,5);
        CP5[:,:,1]=mpoints*B; #bezier control points
        push!(sup,order_plot(CP5,5));=#

    end

    #Corner patches

    for i in 1:length(P[3])
        N=P[3][i][1];
        edge=P[3][i][2]; #face hedge going out from the EV (when present)
        if N==4
            B=get(get(G1S[:masks],"CP",Dict{Int,Any}()),N,acc3bigmatrix(N,"CORNER"));
            mpoints=hm.points[:,P[3][i][3:end]];
            CP3[:,:,1]=mpoints*B;
            order_bezier_cp(gs5,C5,mpoints*B,edge,3);
            #push!(sup,order_plot(CP3,3));
        else
            #if N in collect(keys(get(G1S[:masks],string("CP",S),0)))
            B=get(get(G1S[:masks],string("CP",S),Dict{Int,Any}()),N,G1matrix(N,S,"EVREGBORDER"));
            #else
            #    B=G1matrix(N,S,"EVREGBORDER");
            #end
            mpoints=hm.points[:,P[3][i][3:end]]; #Mesh points corresponding to the neighborhood
            order_bezier_cp(gs5,C5,mpoints*B,edge,5);
            CP5[:,:,1]=mpoints*B; #bezier control points
            #push!(sup,order_plot(CP5,5));
        end
    end

    #Extraordinary boundary patches

    for i in 1:length(P[4])
        face_edge=P[4][i][2];
        l=ccw_edges(hm,face_edge);
        k=P[4][i][1];
        N=2*k;
        if k in collect(keys(get(G1S[:masks],"BOUNDARY",0)))
            B=get(get(G1S[:masks],"BOUNDARY",0),k,0);
        else
            B=G1matrix(N,S,"BORDEREV");
        end
        B1=degree_elevate_mask(N,"BORDEREVR",k); #3 #right patch
        B2=degree_elevate_mask(N,"BORDEREVL",k); #4 #left patch
        B3=degree_elevate_mask(N,"BORDERFACEEV",k) #5; #middle patches
        mpoints=hm.points[:,P[4][i][3:end]];

        B1=ordermatrix(B1,N,"BEV");
        B1[:,4]=B[:,1];
        B1[:,5]=B[:,2*k-1];
        B1[:,7]=B[:,k];
        B1[:,8]=B[:,3*k-1];
        B1[:,33]=B[:,5*k-3];
        order_bezier_cp(gs5,C5,mpoints*B1,l[1],5);
        CP5[:,:,1]=mpoints*B1; #bezier control points

        #push!(sup,order_plot(CP5,5));

        B2=ordermatrix(B2,N,"BEV");
        B2[:,2]=B[:,k-1];
        B2[:,3]=B[:,2*k-2];
        B2[:,5]=B[:,3*k-2];
        B2[:,6]=B[:,5*k-4];
        B2[:,17]=B[:,7*k-6];
        order_bezier_cp(gs5,C5,mpoints*B2,l[end],5);
        CP5[:,:,1]=mpoints*B2; #bezier control points
        #push!(sup,order_plot(CP5,5));

        B3=ordermatrix(B3,N,"BEV");
        for i in 2:k-1
            B3[:,2]=B[:,i-1];
            B3[:,4]=B[:,i];
            B3[:,3]=B[:,k+i-2];
            B3[:,7]=B[:,k+i-1];
            B3[:,5]=B[:,2*k+i-2];
            B3[:,6]=B[:,4*k+i-4];
            B3[:,8]=B[:,3*k+i-2];
            B3[:,33]=B[:,5*k+i-4];
            B3[:,17]=B[:,6*k+i-6];
            order_bezier_cp(gs5,C5,mpoints*B3,l[i],5);
            CP5[:,:,1]=mpoints*B3; #bezier control points
            #push!(sup,order_plot(CP5,5));

            for i in 2:36
                B3[:,i]=completeshift(B3[:,i],k,1,"BORDER");
            end
        end
    end

    gs5.ctrpoints=C5;

    return gs5;

end

function g1surface(m::Mesh, S::String="CS-S")
    hm = hmesh(m);
    divideEV(hm);
    g1surface(hm,S);
end

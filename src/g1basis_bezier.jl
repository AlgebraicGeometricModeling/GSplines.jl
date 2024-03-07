export g1basis_bezier


function g1basis_bezier(m::HMesh)

#m=gs.mesh;
l=ccw_edges(m); #list of vertices with their edges
#basis=spzeros(36*nbf(m),g1dimension(m)[1]);
knts = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
local gsbas = GSpline(0, knts, m )
global globcol=1;

bigrow=Float64[];
bigcol=Float64[];
bigval=Float64[];

for i in 1:length(l) #I go through all the vertices
    e=l[i];
    N=length(e); #number of edges around and vertex valence

    if N==1 #corner
        g1basiscp=g1basis_Corners();
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            #val=findnz(g1basiscp[:,k])[2];
            j=div(row[1],36)+1
            v=e[j];
            h=div(row[1]%36,6)+1;
            w=row[1]%6;
            #basis[idx(gsbas,v,w,h),globcol]=val[1]

            push!(bigrow,idx(gsbas,v,w,h));
            push!(bigcol,globcol);
            push!(bigval,val[1]);
            globcol+=1
        end
    elseif N!=4 && opp(m,e[1])!=0 #inner EV  e[j] is the EV hedge

        g1basiscp=g1basis_EV(N);
        col=size(g1basiscp)[2];
        for k in 1:col
                data=findnz(g1basiscp[:,k]);
                row=data[1];
                val=data[2];
                for i in 1:length(row)
                j=div(row[i],36)+1
                v=e[j];
                h=div(row[i]%36,6)+1;
                w=row[i]%6;
                #basis[idx(gsbas,v,w,h),globcol]=val[i]

                push!(bigrow,idx(gsbas,v,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[i]);
                end
                globcol+=1
        end

    elseif N==2 && opp(m,e[1])==0 #regular border vertex
        if vertexvalence(m,opp(m,e[2]))[1]!=4
            val=vertexvalence(m,opp(m,e[2]))[1]
            g1basiscp=g1basis_borderRV(val);
        else
            g1basiscp=g1basis_borderRV();
        end
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            for i in 1:length(row)
                j=div(row[i],36)+1
                v=e[j];
                h=div(row[i]%36,6)+1;
                w=row[i]%6;
                #basis[idx(gsbas,v,w,h),globcol]=val[i]

                push!(bigrow,idx(gsbas,v,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[i]);
            end
            globcol+=1
        end
    elseif N==4 && opp(m,e[1])!=0 #inner regular vertex
        check=[];
        for s in 1:4
            L=ccw_edges(m,next(m,e[s]))
            if opp(m,L[1])==0 && length(L)>2 #border EV
                push!(check,length(L));
            elseif opp(m,L[1])!=0 && length(L)!=4 #inner EV
                push!(check,length(L));
            else
                push!(check,0);
            end
        end
        g1basiscp=g1basis_RV(check);
        col=size(g1basiscp)[2];
        for k in 1:col
                data=findnz(g1basiscp[:,k]);
                row=data[1];
                val=data[2];
                for i in 1:length(row)
                j=div(row[i],36)+1
                v=e[j];
                h=div(row[i]%36,6)+1;
                w=row[i]%6;
                #basis[idx(gsbas,v,w,h),globcol]=val[i]

                push!(bigrow,idx(gsbas,v,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[i]);
                end
                globcol+=1
        end
    end
end


for i in 1:nbe(m) #I go through all the half edges
    if opp(m,i)!=0 && is_singular_edge(m,i) #inner singular hedge
        g1basiscp=g1basis_EVedge();
        data=findnz(g1basiscp[:,1]);
        row=data[1];
        val=data[2];
        h=div(row[1]%36,6)+1;
        w=row[1]%6;
        #basis[idx(gsbas,i,w,h),globcol]=val[1]

        push!(bigrow,idx(gsbas,i,w,h));
        push!(bigcol,globcol);
        push!(bigval,val[1]);

        h=div(row[2]%36,6)+1;
        w=row[2]%6;
        #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[2]

        push!(bigrow,idx(gsbas,next(m,opp(m,i)),w,h));
        push!(bigcol,globcol);
        push!(bigval,val[2]);

        globcol+=1
    elseif opp(m,i)!=0 && !is_singular_edge(m,i) #inner regular hedge
        g1basiscp=g1basis_Redge();
        col=size(g1basiscp)[2];
        data=findnz(g1basiscp[:,1]);
        row=data[1];
        val=data[2];
        h=div(row[1]%36,6)+1;
        w=row[1]%6;
        #basis[idx(gsbas,i,w,h),globcol]=val[1]

        push!(bigrow,idx(gsbas,i,w,h));
        push!(bigcol,globcol);
        push!(bigval,val[1]);

        h=div(row[2]%36,6)+1;
        w=row[2]%6;
        #basis[idx(gsbas,i,w,h),globcol]=val[2]

        push!(bigrow,idx(gsbas,i,w,h));
        push!(bigcol,globcol);
        push!(bigval,val[2]);

        h=div(row[3]%36,6)+1;
        w=row[3]%6;
        #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[3]

        push!(bigrow,idx(gsbas,next(m,opp(m,i)),w,h));
        push!(bigcol,globcol);
        push!(bigval,val[3]);

        globcol+=1
        data=findnz(g1basiscp[:,2]);
        row=data[1];
        val=data[2];
        h=div(row[1]%36,6)+1;
        w=row[1]%6;
        #basis[idx(gsbas,i,w,h),globcol]=val[1]

        push!(bigrow,idx(gsbas,i,w,h));
        push!(bigcol,globcol);
        push!(bigval,val[1]);

        h=div(row[2]%36,6)+1;
        w=row[2]%6;
        #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[2]

        push!(bigrow,idx(gsbas,next(m,opp(m,i)),w,h));
        push!(bigcol,globcol);
        push!(bigval,val[2]);

        h=div(row[3]%36,6)+1;
        w=row[3]%6;
        #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[3]

        push!(bigrow,idx(gsbas,next(m,opp(m,i)),w,h));
        push!(bigcol,globcol);
        push!(bigval,val[3]);

        globcol+=1
    elseif opp(m,i)==0 && !is_singular_edge(m,i) #border regular hedge
        g1basiscp=g1basis_RBedge();
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            h=div(row[1]%36,6)+1;
            w=row[1]%6;
            #basis[idx(gsbas,i,w,h),globcol]=val[1]

            push!(bigrow,idx(gsbas,i,w,h));
            push!(bigcol,globcol);
            push!(bigval,val[1]);

            globcol+=1
        end
    end
end

for i in 1:nbf(m)
    g1basiscp=g1basis_Faces();
    col=size(g1basiscp)[2];
    for k in 1:col
        data=findnz(g1basiscp[:,k]);
        row=data[1];
        val=data[2];
        h=div(row[1]%36,6)+1;
        w=row[1]%6;
        #basis[idx(gsbas,hm.faces[i],w,h),globcol]=val[1]

        push!(bigrow,idx(gsbas,m.faces[i],w,h));
        push!(bigcol,globcol);
        push!(bigval,val[1]);

        globcol+=1
    end
end

return basis=sparse(bigrow,bigcol,bigval)



end

function g1basis_bezier(m::Mesh)
    hm=hmesh(m);
    g1basis_bezier(hm);
end

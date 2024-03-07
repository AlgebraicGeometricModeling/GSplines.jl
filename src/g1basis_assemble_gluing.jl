function assemble_g1basis_gluing(hm::HMesh, knts::Array, subd_step::Int64)

l=ccw_edges(hm); #list of vertices with their edges
#basis=spzeros(36*nbf(m),g1dimension(m)[1]);
#knts = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5 ,1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
m, d = dim_deg(knts)
local gsbas = GSpline(0, knts, hm)
global globcol=1;

bigrow=Float64[];
bigcol=Float64[];
bigval=Float64[];

for i in 1:length(l) #I go through all the vertices
    e=l[i];
    N=length(e); #number of edges around and vertex valence

    if N==1 #corner
        g1basiscp=g1basis_Corners_kn(knts);
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            #val=findnz(g1basiscp[:,k])[2];
            j=div(row[1],m^2)+1
            v=e[j];
            h=div(row[1]%m^2,m)+1;
            w=row[1]%m;
            #basis[idx(gsbas,v,w,h),globcol]=val[1]

            push!(bigrow,idx(gsbas,v,w,h));
            push!(bigcol,globcol);
            push!(bigval,val[1]);
            globcol+=1
        end
    elseif N!=4 && opp(hm,e[1])!=0 #inner EV  e[j] is the EV hedge

        glu0=[2*cos(2*pi/N),0,0];
        glu=glob_gluing_subd(glu0,subd_step);

        #g1basiscp=g1basis_EV_gluing(knts,glu[1],N)

        g1basiscp=g1basis_EV(N);
        g1basiscp=bezier_to_spline(g1basiscp,knts,N,"innervertex");

        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            for i in 1:length(row)
                j=div(row[i],m^2)+1
                #println(j)
                v=e[j];
                h=div(row[i]%m^2,m)+1;
                w=row[i]%m;
                #basis[idx(gsbas,v,w,h),globcol]=val[i]
                push!(bigrow,idx(gsbas,v,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[i]);
            end
            globcol+=1
        end

        if subd_step>0
            m0=6;

            for j in 1:2^subd_step -1

                gluL=glu[j];
                #println(gluL)
                gluR=glu[j+1];
                #println(gluR)
                p=m0+(m0-1)*(j-1);

                bas=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]
                for b in 1:length(bas)
                g1basiscp=g1basis_noncrossingRV(knts, gluL, gluR,p,1,N)[:,bas[b]];
                #g1basiscp=g1basis_noncrossingRV_alternative(knts, gluL, gluR,p,1,N);

                col=1#size(g1basiscp)[2];
                for k in 1:col
                    data=findnz(g1basiscp[:,k]);
                    row=data[1];
                    val=data[2];
                    for i in 1:length(row)
                        j=div(row[i],m^2)+1
                        #println(j)
                        v=e[j];
                        h=div(row[i]%m^2,m)+1;
                        w=row[i]%m;
                        #basis[idx(gsbas,v,w,h),globcol]=val[i]
                        push!(bigrow,idx(gsbas,v,w,h));
                        push!(bigcol,globcol);
                        push!(bigval,val[i]);
                    end

                    globcol+=1
                end
                end

            end

        end


    elseif N==2 && opp(hm,e[1])==0 #regular border vertex
            #=#if vertexvalence(hm,opp(hm,e[2]))[1]!=4
                val=vertexvalence(hm,opp(hm,e[2]))[1];

                glu0=[2*cos(2*pi/val),0,0];

                glu=glob_gluing_subd(glu0,subd_step);

                g1basiscp=g1basis_borderRV_gluing(knts,glu[end]);=#

                if vertexvalence(hm,opp(hm,e[2]))[1]!=4
                    val=vertexvalence(hm,opp(hm,e[2]))[1]
                    #g1basiscp=g1basis_RBV_kn(knts,val);

                    g1basiscp=g1basis_borderRV(val);
                    g1basiscp=bezier_to_spline(g1basiscp,knts,N,"bordervertex");
                else
                    #g1basiscp=g1basis_RBV_kn(knts);

                    g1basiscp=g1basis_borderRV();
                    g1basiscp=bezier_to_spline(g1basiscp,knts,N,"bordervertex");
                end

        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            for i in 1:length(row)
                j=div(row[i],m^2)+1
                v=e[j];
                h=div(row[i]%m^2,m)+1;
                w=row[i]%m;
                #basis[idx(gsbas,v,w,h),globcol]=val[i]

                push!(bigrow,idx(gsbas,v,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[i]);
            end
            globcol+=1
        end
    elseif N==4 && opp(hm,e[1])!=0 #inner regular vertex
        check=[];
        nbrv=0;

        glu0=[2*cos(2*pi/N),0,0];
        glu=glob_gluing_subd(glu0,subd_step);

        g1basiscp=g1basis_EV_gluing(knts,glu[1],N)

        for s in 1:4
            L=ccw_edges(hm,next(hm,e[s]))
            if opp(hm,L[1])==0 && length(L)>2 #border EV
                push!(check,length(L));
            elseif opp(hm,L[1])!=0 && length(L)!=4 #inner EV
                push!(check,length(L));
            else
                push!(check,0);
                nbrv+=1;
            end
        end

        if nbrv==4 #no EVs around, no knot insertion needed
            g1basiscp=g1basis_RVs_coeff_kn(knts)
            #g1basiscp=g1basis_RV(check);
        else
            g1basiscp=g1basis_RV_kn(knts,check)

            #g1basiscp=g1basis_RV(check);
            #g1basiscp=bezier_to_spline(g1basiscp,knts,N,"innervertex");
        end

        col=size(g1basiscp)[2];

        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            for i in 1:length(row)
                j=div(row[i],m^2)+1
                v=e[j];
                h=div(row[i]%m^2,m)+1;
                w=row[i]%m;
                #basis[idx(gsbas,v,w,h),globcol]=val[i]
                push!(bigrow,idx(gsbas,v,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[i]);
            end
            globcol+=1
        end

        if subd_step>0
            m0=6;

            for j in 1:2^subd_step -1

                gluL=glu[j];
                #println(gluL)
                gluR=glu[j+1];
                #println(gluR)
                p=m0+(m0-1)*(j-1);

                g1basiscp=g1basis_noncrossingRV_alternative(knts, gluL, gluR,p,1,N);

                col=size(g1basiscp)[2];
                for k in 1:col
                    data=findnz(g1basiscp[:,k]);
                    row=data[1];
                    val=data[2];
                    for i in 1:length(row)
                        j=div(row[i],m^2)+1
                        #println(j)
                        v=e[j];
                        h=div(row[i]%m^2,m)+1;
                        w=row[i]%m;
                        #basis[idx(gsbas,v,w,h),globcol]=val[i]
                        push!(bigrow,idx(gsbas,v,w,h));
                        push!(bigcol,globcol);
                        push!(bigval,val[i]);
                    end

                    globcol+=1
                end

            end
        end
    end
end

list=[]

for i in 1:nbe(hm) #I go through all the half edges
    if opp(hm,i) in list
        continue
    end
        if opp(hm,i)!=0 && is_singular_edge(hm,i) #inner singular hedge

            #g1basiscp=g1basis_EVedge_gluing(knts,subd_step);
            g1basiscp=g1basis_EVedge_kn(knts);
            col=size(g1basiscp)[2];
            for k in 1:col
                data=findnz(g1basiscp[:,k]);
                row=data[1];
                val=data[2];
                h=div(row[1]%m^2,m)+1;
                w=row[1]%m;
                #basis[idx(gsbas,i,w,h),globcol]=val[1]
                push!(bigrow,idx(gsbas,i,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[1]);

                h=div(row[2]%m^2,m)+1;
                w=row[2]%m;
                #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[2]

                push!(bigrow,idx(gsbas,next(hm,opp(hm,i)),w,h));
                push!(bigcol,globcol);
                push!(bigval,val[2]);

                globcol+=1
            end

            push!(list,i)

    elseif opp(hm,i)!=0 && !is_singular_edge(hm,i) #inner regular hedge

        g1basiscp=g1basis_Redge_gluing(knts,subd_step);
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            if k%2!=0
                row=data[1];
                val=data[2];
                h=div(row[1]%m^2,m)+1;
                w=row[1]%m;
                #basis[idx(gsbas,i,w,h),globcol]=val[1]

                push!(bigrow,idx(gsbas,i,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[1]);

                h=div(row[2]%m^2,m)+1;
                w=row[2]%m;
                #basis[idx(gsbas,i,w,h),globcol]=val[2]

                push!(bigrow,idx(gsbas,i,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[2]);

                h=div(row[3]%m^2,m)+1;
                w=row[3]%m;
                #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[3]

                push!(bigrow,idx(gsbas,next(hm,opp(hm,i)),w,h));
                push!(bigcol,globcol);
                push!(bigval,val[3]);

                globcol+=1
            else
                row=data[1];
                val=data[2];
                h=div(row[1]%m^2,m)+1;
                w=row[1]%m;
                #basis[idx(gsbas,i,w,h),globcol]=val[1]

                push!(bigrow,idx(gsbas,i,w,h));
                push!(bigcol,globcol);
                push!(bigval,val[1]);

                h=div(row[2]%m^2,m)+1;
                w=row[2]%m;
                #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[2]

                push!(bigrow,idx(gsbas,next(hm,opp(hm,i)),w,h));
                push!(bigcol,globcol);
                push!(bigval,val[2]);

                h=div(row[3]%m^2,m)+1;
                w=row[3]%m;
                #basis[idx(gsbas,next(m,opp(m,i)),w,h),globcol]=val[3]

                push!(bigrow,idx(gsbas,next(hm,opp(hm,i)),w,h));
                push!(bigcol,globcol);
                push!(bigval,val[3]);

                globcol+=1
            end
        end
        push!(list,i)
    elseif opp(hm,i)==0 && !is_singular_edge(hm,i) #border regular hedge
        g1basiscp=g1basis_RBedge_kn(knts);
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            h=div(row[1]%m^2,m)+1;
            w=row[1]%m;
            #basis[idx(gsbas,i,w,h),globcol]=val[1]

            push!(bigrow,idx(gsbas,i,w,h));
            push!(bigcol,globcol);
            push!(bigval,val[1]);

            globcol+=1
        end
        push!(list,i)
    end
end


    for i in 1:nbf(hm)
        g1basiscp=g1basis_Faces_kn(knts);
        col=size(g1basiscp)[2];
        for k in 1:col
            data=findnz(g1basiscp[:,k]);
            row=data[1];
            val=data[2];
            h=div(row[1]%m^2,m)+1;
            w=row[1]%m;
            #basis[idx(gsbas,hm.faces[i],w,h),globcol]=val[1]

            push!(bigrow,idx(gsbas,hm.faces[i],w,h));
            push!(bigcol,globcol);
            push!(bigval,val[1]);

            globcol+=1
        end
    end
    #println(bigrow)
return basis=sparse(bigrow,bigcol,bigval);

end

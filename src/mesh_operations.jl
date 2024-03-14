export divideEV

function kindofface(m::HMesh)
F=nbf(m);
innerf=Int64[];
borderf=Int64[];
cornerf=Int64[];
for i in 1:F
    e=Int64[];
    z=0;
    push!(e,i);
    fe=edges_on_face(m,i);
    for j in 1:convert(Int64,length(fe))
        if opp(m,fe[j])==0 #border edge
            push!(e,0);
            z=z+1;
        else
            push!(e,1);
        end
    end

    if z==0 #interior face
        push!(innerf,i);

    elseif z==1 #border face
        push!(borderf,i);

    else #corner face
        push!(cornerf,i);

    end

end

    #Isolate border patches with EV

    EVborderf=Int64[];
    REGborderf=Int64[];
    for i in 1:length(borderf)
        e=m.faces[borderf[i]]; #face hedge
        while opp(m,e)!=0 #select border hedge with no opposite
            e=next(m,e);
        end
        if opp(m,prev(m,opp(m,prev(m,e))))==0 && opp(m,next(m,opp(m,next(m,e))))==0  #case of a regular border patch with k=2
        push!(REGborderf,borderf[i]);
        else
        push!(EVborderf,e);
        end
    end

    #Isolate inner patches with a border EV vertex

    EVinnerf=Int64[];
    REGinnerf=Int64[];
    for j in 1:length(innerf)
        z=0;
        e=m.faces[innerf[j]]; #face hedge
        for i in 1:4
            f=ccw_edges(m,e);
            if opp(m,prev(m,f[end]))!=0
                e=next(m,e);
            else
                push!(EVinnerf,e);
                z=z+1;
            end
        end
        if z==0
            push!(REGinnerf,innerf[j]);
        end
    end


    return REGinnerf,REGborderf,cornerf,EVinnerf,EVborderf

end


#--------------------------------------------------------------------------------------

function vertexvalence(m::HMesh, e::Int64) #find the edge corresponding to an ev in a REGULAR INNER face
    #=if b==0
        l=ccw_edges(m,e)
        if opp(m,prev(m,l[end]))==0 #border EV
            while opp(m,e)!=0
                e=next(m,opp(m,e))
            end
            v=length(ccw_edges(m,e))
            return v,e
    else=#
    ef=e;
    for i in 1:4
      ei=e;
      c=1;
      if opp(m,e)!=0 && opp(m,next(m,opp(m,e)))!=0
          while ei != next(m,opp(m,e))
              e=next(m,opp(m,e));
              c=c+1; #vertex valence
          end
          if c != 4
              return c,ei
          end
      end
              e=next(m,ei);
      end
      return 4,ef
end

#-------------------------------------------------------------------------------------


function neigh(m::HMesh)

    inner=kindofface(m)[1];
    border=kindofface(m)[2];
    corner=kindofface(m)[3];
    inner_points=Int64[];
    border_points=Int64[];
    corner_points=Int64[];

    for i in 1:length(inner) #inner faces, complete ones
        fneigh=Int64[];
        sneigh=Int64[];
        Aneigh=Int64[];
        Bneigh=Int64[];
        Cneigh=Int64[];
        Dneigh=Int64[];
        innerneighpoints=Int64[];
        e=m.faces[inner[i]]; #face hedge
        V=vertexvalence(m,e);
        if V[1] !=4
            e=V[2];
        end

        for i in 1:4
            fneigh=Int64[];
            sneigh=Int64[];
            Aneigh=Int64[];
            Bneigh=Int64[];
            Cneigh=Int64[];
            Dneigh=Int64[];
        push!(fneigh,V[1]); #vertex valence
        push!(fneigh,e);
        push!(fneigh,e); #center mask
        ne=next(m,e);
            try
        push!(fneigh,ne); #right 1st neigh point
        push!(sneigh,next(m,ne)); #1st value of 2nd neigh
        push!(Aneigh,next(m,next(m,opp(m,ne)))); #1st value of Aneigh
        push!(Bneigh,next(m,next(m,next(m,opp(m,ne))))); #1st value of Bneigh
        push!(Cneigh,next(m,next(m,opp(m,next(m,ne))))); #1st value of Cneigh
        push!(Dneigh,prev(m,opp(m,next(m,opp(m,next(m,ne)))))); #1st value of Dneigh

        edge1=ne;

        while edge1 != next(m,opp(m,prev(m,(prev(m,ne)))))

            ne=next(m,opp(m,prev(m,(prev(m,ne)))));
            push!(fneigh,ne);
            push!(sneigh,next(m,ne));
            push!(Aneigh,next(m,next(m,opp(m,ne))));
            push!(Bneigh,next(m,next(m,next(m,opp(m,ne)))));
            push!(Cneigh,next(m,next(m,opp(m,next(m,ne)))));
            push!(Dneigh,prev(m,opp(m,next(m,opp(m,next(m,ne))))));
        end
            catch ne
                if isa(ne,BoundsError)
                    e=next(m,e);
                end
            end
        end


        innerneighpoints=hcat(fneigh',sneigh',Aneigh',Bneigh',Dneigh',Cneigh');

        for i in 3:convert(Int64,length(innerneighpoints))
            innerneighpoints[i]=ptidx_of(m,innerneighpoints[i]);
        end
        inner_points=vcat(inner_points,[innerneighpoints]);
    end

    for i in 1:length(border) #border faces

        fneigh=Int64[];
        sneigh=Int64[];
        borderneighpoints=Int64[];
        e=m.faces[border[i]];
        V=vertexvalence(m,e);
        if V[1] !=4
            val=V[1];
            e=V[2];
        else
            while opp(m,next(m,next(m,e)))!=0
                e=next(m,e);
            end
            val=4;
        end

        push!(fneigh,val); #vertex valence
        push!(fneigh,e); #center mask
        push!(fneigh,e);
        ne=next(m,e);
        push!(fneigh,ne); #right 1st neigh point
        push!(sneigh,next(m,ne)); #1st value of 2nd neigh

        edge1=ne;

        while edge1 != next(m,opp(m,prev(m,(prev(m,ne)))))

            ne=next(m,opp(m,prev(m,(prev(m,ne)))));
            push!(fneigh,ne);
            push!(sneigh,next(m,ne));

        end

        if opp(m,next(m,next(m,e)))==0
            q0=next(m,next(m,opp(m,fneigh[4])));
            q0b=next(m,q0);
            q3b=prev(m,opp(m,prev(m,q0)));
        elseif opp(m,next(m,e))==0
            q0=prev(m,prev(m,opp(m,prev(m,e))));
            q0b=prev(m,q0);
            q3b=next(m,next(m,opp(m,e)));
        end

        borderneighpoints=hcat(fneigh',sneigh',[q0,q0b,q3b]');

        for i in 3:convert(Int64,length(borderneighpoints))
            borderneighpoints[i]=ptidx_of(m,borderneighpoints[i]);
        end
        border_points=vcat(border_points,[borderneighpoints]);
    end


    for i in 1:length(corner); #border faces, no problems with 1st and 2nd neighborhood

        fneigh=Int64[];
        sneigh=Int64[];
        cornerneighpoints=Int64[];
        e=m.faces[corner[i]];
        V=vertexvalence(m,e);
        if V[1] !=4
            val=V[1];
            e=V[2];
        else
            while opp(m,e)==0 || opp(m,next(m,e))!=0
                e=next(m,e);
            end
            val=4;
        end


        push!(fneigh,val); #vertex valence
        push!(fneigh,e); #center mask
        push!(fneigh,e);
        ne=next(m,e);
        push!(fneigh,ne); #right 1st neigh point
        push!(sneigh,next(m,ne)); #1st value of 2nd neigh

        edge1=ne;

        while edge1 != next(m,opp(m,prev(m,(prev(m,ne)))))

            ne=next(m,opp(m,prev(m,(prev(m,ne)))));
            push!(fneigh,ne);
            push!(sneigh,next(m,ne));

        end

        cornerneighpoints=hcat(fneigh',sneigh');

        for i in 3:convert(Int64,length(cornerneighpoints))
            cornerneighpoints[i]=ptidx_of(m,cornerneighpoints[i]);
        end
        corner_points=vcat(corner_points,[cornerneighpoints]);
    end


    #Neighborhood for borders EVs

    EVborder=kindofface(m)[5];
    rightEVborder=Int64[];
    for i in 1:length(EVborder)
        if length(ccw_edges(m,EVborder[i]))!=2
            push!(rightEVborder,EVborder[i])
        end
    end

    fakeneigh_points=Int64[];

        for i in 1:length(rightEVborder)
            fneigh=Int64[];
            sneigh=Int64[];
            Aneigh=Int64[];
            Bneigh=Int64[];
            Cneigh=Int64[];
            Dneigh=Int64[];
            boundEVneighpoints=Int64[];
            e=rightEVborder[i];

            k=length(ccw_edges(m,e));
            push!(fneigh,k); #boundary patches valance
            push!(fneigh,e); #center mask
            push!(fneigh,e); #lo metto due volte per ricordarmi l'hedge e non convertirlo in punto per dopo
            ne=next(m,e);
            push!(fneigh,ne); #right 1st neigh point
            push!(sneigh,next(m,ne)); #1st value of 2nd neigh
            push!(Aneigh,next(m,next(m,opp(m,ne)))); #1st value of Aneigh
            push!(Bneigh,next(m,next(m,next(m,opp(m,ne))))); #1st value of Bneigh
            push!(Cneigh,next(m,next(m,opp(m,next(m,ne))))); #1st value of Cneigh
            push!(Dneigh,prev(m,opp(m,next(m,opp(m,next(m,ne)))))); #1st value of Dneigh
            for i in 1:k-1
                ne=next(m,opp(m,prev(m,(prev(m,ne)))));
                push!(fneigh,ne);
                push!(sneigh,next(m,ne));
                push!(Aneigh,next(m,next(m,opp(m,ne))));
                push!(Bneigh,next(m,next(m,next(m,opp(m,ne)))));
                push!(Cneigh,next(m,next(m,opp(m,next(m,ne)))));
                push!(Dneigh,prev(m,opp(m,next(m,opp(m,next(m,ne))))));
            end
            push!(fneigh,prev(m,prev(m,ne))); #I complete the neighborhoods
            push!(Aneigh,prev(m,opp(m,prev(m,prev(m,prev(m,ne))))));

            #boundEVneighpoints=hcat(fneigh',fneigh[5:end-1]',sneigh',sneigh',Aneigh',Aneigh[2:end-1]',Bneigh',Bneigh',Dneigh',Dneigh',Cneigh',Cneigh');
            boundEVneighpoints=hcat(fneigh',sneigh',Aneigh',Bneigh',Dneigh',Cneigh');

            for i in 3:convert(Int64,length(boundEVneighpoints))
                boundEVneighpoints[i]=ptidx_of(m,boundEVneighpoints[i]);
            end
            fakeneigh_points=vcat(fakeneigh_points,[boundEVneighpoints]);
        end

        return inner_points,border_points,corner_points,fakeneigh_points
end





function neighacc3(m::HMesh)

    inner=kindofface(m)[1];
    border=kindofface(m)[2];
    corner=kindofface(m)[3];
    inner_points=Int64[];
    border_points=Int64[];
    corner_points=Int64[];

    for i in 1:length(inner) #inner faces, complete ones
        fneigh=Int64[];
        sneigh=Int64[];

        innerneighpoints=Int64[];
        e=m.faces[inner[i]]; #face hedge
        V=vertexvalence(m,e);
        if V[1] !=4
            e=V[2];
        end

        for i in 1:4
            fneigh=Int64[];
            sneigh=Int64[];

        push!(fneigh,V[1]); #vertex valence
        push!(fneigh,e);
        push!(fneigh,e); #center mask
        ne=next(m,e);
            try
        push!(fneigh,ne); #right 1st neigh point
        push!(sneigh,next(m,ne)); #1st value of 2nd neigh

        edge1=ne;

        while edge1 != next(m,opp(m,prev(m,(prev(m,ne)))))

            ne=next(m,opp(m,prev(m,(prev(m,ne)))));
            push!(fneigh,ne);
            push!(sneigh,next(m,ne));

        end
            catch ne
                if isa(ne,BoundsError)
                    e=next(m,e);
                end
            end
        end


        innerneighpoints=hcat(fneigh',sneigh');

        for i in 3:convert(Int64,length(innerneighpoints))
            innerneighpoints[i]=ptidx_of(m,innerneighpoints[i]);
        end
        inner_points=vcat(inner_points,[innerneighpoints]);
    end

    for i in 1:length(border) #border faces

        fneigh=Int64[];
        sneigh=Int64[];
        borderneighpoints=Int64[];
        e=m.faces[border[i]];
        V=vertexvalence(m,e);
        if V[1] !=4
            val=V[1];
            e=V[2];
        else
            while opp(m,next(m,next(m,e)))!=0
                e=next(m,e);
            end
            val=4;
        end

        push!(fneigh,val); #vertex valence
        push!(fneigh,e); #center mask
        push!(fneigh,e);
        ne=next(m,e);
        push!(fneigh,ne); #right 1st neigh point
        push!(sneigh,next(m,ne)); #1st value of 2nd neigh

        edge1=ne;

        while edge1 != next(m,opp(m,prev(m,(prev(m,ne)))))

            ne=next(m,opp(m,prev(m,(prev(m,ne)))));
            push!(fneigh,ne);
            push!(sneigh,next(m,ne));

        end


        borderneighpoints=hcat(fneigh',sneigh');

        for i in 3:convert(Int64,length(borderneighpoints))
            borderneighpoints[i]=ptidx_of(m,borderneighpoints[i]);
        end
        border_points=vcat(border_points,[borderneighpoints]);
    end


    for i in 1:length(corner); #border faces, no problems with 1st and 2nd neighborhood

        fneigh=Int64[];
        sneigh=Int64[];
        cornerneighpoints=Int64[];
        e=m.faces[corner[i]];
        V=vertexvalence(m,e);
        if V[1] !=4
            val=V[1];
            e=V[2];
        else
            while opp(m,e)==0 || opp(m,next(m,e))!=0
                e=next(m,e);
            end
            val=4;
        end


        push!(fneigh,val); #vertex valence
        push!(fneigh,e); #center mask
        push!(fneigh,e);
        ne=next(m,e);
        push!(fneigh,ne); #right 1st neigh point
        push!(sneigh,next(m,ne)); #1st value of 2nd neigh

        edge1=ne;

        while edge1 != next(m,opp(m,prev(m,(prev(m,ne)))))

            ne=next(m,opp(m,prev(m,(prev(m,ne)))));
            push!(fneigh,ne);
            push!(sneigh,next(m,ne));

        end

        cornerneighpoints=hcat(fneigh',sneigh');

        for i in 3:convert(Int64,length(cornerneighpoints))
            cornerneighpoints[i]=ptidx_of(m,cornerneighpoints[i]);
        end
        corner_points=vcat(corner_points,[cornerneighpoints]);
    end


    #Neighborhood for borders EVs

    EVborder=kindofface(m)[5];
    rightEVborder=Int64[];
    for i in 1:length(EVborder)
        if length(ccw_edges(m,EVborder[i]))!=2
            push!(rightEVborder,EVborder[i])
        end
    end

    fakeneigh_points=Int64[];

        for i in 1:length(rightEVborder)
            fneigh=Int64[];
            sneigh=Int64[];
            Aneigh=Int64[];
            Bneigh=Int64[];
            Cneigh=Int64[];
            Dneigh=Int64[];
            boundEVneighpoints=Int64[];
            e=rightEVborder[i];

            k=length(ccw_edges(m,e));
            push!(fneigh,k); #boundary patches valance
            push!(fneigh,e); #center mask
            push!(fneigh,e); #lo metto due volte per ricordarmi l'hedge e non convertirlo in punto per dopo
            ne=next(m,e);
            push!(fneigh,ne); #right 1st neigh point
            push!(sneigh,next(m,ne)); #1st value of 2nd neigh
            push!(Aneigh,next(m,next(m,opp(m,ne)))); #1st value of Aneigh
            push!(Bneigh,next(m,next(m,next(m,opp(m,ne))))); #1st value of Bneigh
            push!(Cneigh,next(m,next(m,opp(m,next(m,ne))))); #1st value of Cneigh
            push!(Dneigh,prev(m,opp(m,next(m,opp(m,next(m,ne)))))); #1st value of Dneigh
            for i in 1:k-1
                ne=next(m,opp(m,prev(m,(prev(m,ne)))));
                push!(fneigh,ne);
                push!(sneigh,next(m,ne));
                push!(Aneigh,next(m,next(m,opp(m,ne))));
                push!(Bneigh,next(m,next(m,next(m,opp(m,ne)))));
                push!(Cneigh,next(m,next(m,opp(m,next(m,ne)))));
                push!(Dneigh,prev(m,opp(m,next(m,opp(m,next(m,ne))))));
            end
            push!(fneigh,prev(m,prev(m,ne))); #I complete the neighborhoods
            push!(Aneigh,prev(m,opp(m,prev(m,prev(m,prev(m,ne))))));

            #boundEVneighpoints=hcat(fneigh',fneigh[5:end-1]',sneigh',sneigh',Aneigh',Aneigh[2:end-1]',Bneigh',Bneigh',Dneigh',Dneigh',Cneigh',Cneigh');
            boundEVneighpoints=hcat(fneigh',sneigh',Aneigh',Bneigh',Dneigh',Cneigh');

            for i in 3:convert(Int64,length(boundEVneighpoints))
                boundEVneighpoints[i]=ptidx_of(m,boundEVneighpoints[i]);
            end
            fakeneigh_points=vcat(fakeneigh_points,[boundEVneighpoints]);
        end

        return inner_points,border_points,corner_points,fakeneigh_points
end

#----------------------------------------------------------------------------------------
function isEV(m::HMesh,e::Int64)
    k=ccw_edges(m,e);
    if length(k)>=3 && opp(m,k[1])==0 #first right hedge of a border EV
        return 1
    elseif opp(m,prev(m,k[end]))==0 #middle and last border EV hedge
        v=k[end];
        while opp(m,v)!=0
            v=next(m,opp(m,v));
        end
        if length(ccw_edges(m,v))>=3
            return 1
        elseif length(ccw_edges(m,v))==2
            return 2
        else
            return 0
        end
    elseif length(k)==3 || length(k)>=5 #inner EV hedge
            return 1
    else
        return 0
    end
end

function divideEV(m::HMesh)
    k=filter(!isempty,ccw_edges(m));
    for i in 1:length(k)
        vertex=k[i];
        if isEV(m,vertex[1])==1 #it is an EV
            for j in 1:length(vertex)
                p=next(m,vertex[j])
                pp=next(m,p)
                if isEV(m,p)==1 || isEV(m,pp)==1
                    @warn("The mesh has two adjacent EVs. One step of Catmull-Clark subdivision will be performed to continue the computation.")
                    cc_subdivide!(m);
                    return
                #elseif isEV(m,p)==2 && opp(m,p)==0
                #    @warn("The mesh has an EVs linked to a boundary vertex. One step of Catmull-Clark subdivision will be performed to continue the computation.")
                #    cc_subdivide!(m);
                    #return
                end
            end
        end
    end
end

function nbevs(m::HMesh,p::Int64=1)

    P=neigh(m);

    if p==1
        println("Number of vertices: ", nbv(m))
        println("Number of faces: ", nbf(m))
    end

    ev=[];
    evsb=[];

    for i in 1:length(P[1])
        N=P[1][i][1]; #EV valence
        if N!=4
            push!(ev,N);
        end
    end

    for i in 1:length(P[3])
        N=P[3][i][1]; #EV valence
        if N!=4
            push!(ev,N);
        end
    end

    for i in 1:length(P[4])
        kk=P[4][i][1];
        push!(evsb,kk);
    end

    nev=[];
    bevpatches=[];
    c=counter(ev);
    k=collect(keys(c));
    v=collect(values(c));
    bev=counter(evsb);
    nbev=collect(keys(bev));
    vbev=collect(values(bev));

    for i in 1:length(k)
        if p==1
            print("Valence ", k[i]," EVs: ")
        end
        evs=v[i]/k[i];
        push!(nev,evs);
        push!(bevpatches,k[i]*evs);
        if p==1
            println(convert(Int64,evs));
        end
    end

    for i in 1:length(evsb)
        if p==1
            print("Valence ", nbev[i]," BEVs: ")
            println(convert(Int64,vbev[i]));
        end
        push!(bevpatches,nbev[i]*vbev[i]);
    end

    if p==1
        if length(ev)!=0
            S=convert(Int64,sum(nev));
            println("Total Inner EVs: ",S )
        end

        if length(evsb)!=0
            KK=length(evsb);
            println("Total Boundary EVs: ",KK )
        end

        if length(ev)!=0 && length(evsb)!=0
            println("Total Number of EVs: ", S+KK)
        end
    end
    #return sum(bevpatches)
end

function nbevs(m::Mesh,p::Int64=1)
    hm=hmesh(m);
    nbevs(hm,p)
end
# the mesh is assumed FIXED ( fix_hedge )

function edge_of(m::HMesh, v::Int64)
    for e in 1:nbe(m)
        if ptidx_of(m,e)==v # && opp(m,e)==0
            e1 = e;
            if opp(m,e1) == 0
                return e1;
            end
            e1 = next(m,opp(m,e1));
            while opp(m,e1) != 0 && e1 != e
                 e1 = next(m,opp(m,e1));
            end
            return e1;
        end
    end
end

function valence_of(m::HMesh, v::Int64)
    e = edge_of(m,v)
    valence=1;
    e1 = opp(m,prev(m,e));
    if e1 == 0 #CORNER CASE
        return 1;
    end
    while e1 != 0 && e1 != e
         e1 = opp(m,prev(m,e1));
         valence+=1
    end
    return valence;
end


#=function fix_hedge(m::HMesh)
    for i in 1:nbv(m)
        e=edge_of_unsorted(m,i)
        tmp=m.edges[i]
        m.edges[i]=m.edges[e]
        m.edges[e]=tmp
    end
end=#
function is_boundary_edge(m::HMesh,e::Int64)
    return opp(m,e)==0
end

function is_boundary_vertex(m::HMesh,v::Int64)
    e=edge_of(m,v);
    return opp(m,e)==0
end

function nb_boundary_edge(m::HMesh)
    nb=0;
    for i in 1:nbe(m)
        nb+=is_boundary_edge(m,i)
    end
    return nb
end

function nb_edge(m::HMesh)
    return convert(Int64,(nbe(m)+nb_boundary_edge(m))/2)
end

function is_corner(m::HMesh,v::Int64)
    e=edge_of(m,v);
    return opp(m,prev(m,e))==0
end

function is_EV(m::HMesh,v::Int64)
    val=valence_of(m,v);
    if (is_boundary_vertex(m,v))
        return val>2 #border EV
    else
        return val!=4 #inner EV
    end
end

function is_singular_edge(m::HMesh,e::Int64)
    return is_EV(m,m.edges[e].point) || is_EV(m,m.edges[next(m,e)].point)
end


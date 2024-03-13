export redefine_boundary_edges!
"""
This function modifies the half edge data structure of a mesh re-defining some edges to be boundaries. 
    It is done giving in input a list of indices identifying the edges (in the half edge structure) to be set as boundaries.
"""
function redefine_boundary_edges!(hm::HMesh, sharp_edges::Array)
    
list_edges=hm.edges;
nb_edges=length(list_edges);
nb_sh_edges=size(sharp_edges,1);

for i in 1:nb_edges

    v1=list_edges[i].point; #Two vertices of the edge in the vertices list
    v2=ptidx_of(hm,list_edges[i].nxt);

if list_edges[i].opp !=0 #If we are not already considering a boundary edge
    for j in 1:nb_sh_edges

        sh_v1=sharp_edges[j,1];
        sh_v2=sharp_edges[j,2];

        if (v1==sh_v1 && v2==sh_v2) || (v1==sh_v2 && v2==sh_v1)
            opp1=list_edges[i].opp;
            list_edges[i].opp=0;
            list_edges[opp1].opp=0;
        end
    end
end
end

end

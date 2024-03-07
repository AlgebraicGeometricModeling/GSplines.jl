export get_basis_from_sparse
#This function compute the ctrpts in order to plot the G1 basis functions

function get_basis_from_sparse(b::SparseVector,gsurf::GSpline)
    local gsbas = GSpline(0, gsurf.knots, gsurf.mesh)
    row=findnz(b)[1];
    val=findnz(b)[2];
    CP=copy(gsurf.ctrpoints);
    for j in 1:length(row)
        CP[3,row[j]]=val[j];
    end
    gsbas.ctrpoints=CP;
    return gsbas
end
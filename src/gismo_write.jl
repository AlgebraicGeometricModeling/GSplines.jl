export ToGismo

#include("src/index.jl")


## Converts surface to gismo xml format
# gs: GSpline
# fl: filename to write to
function AddSurfaceToGismoFile(gs,g, F::Int64)
    for f in 1:F
        write(g,string("<Geometry type=\"TensorBSpline2\" id=\" ",f-1,"\" >\n"))
        write(g,"<Basis type=\"TensorBSplineBasis2\">\n")
        write(g,"<Basis type=\"BSplineBasis\" index=\" 0\"> \n")
        write(g,"<KnotVector degree=\" ")
        write(g,string(gs.degree))
        write(g,"\"> ")
        #for i in 1:length(gs.knots)            write(g,string(gs.knots[i]," "))        end
        write(g,string(gs.knots," </KnotVector> \n"))
        write(g,"</Basis>\n")
        write(g,"<Basis type=\"BSplineBasis\" index=\" 1\" >\n")
        write(g,"<KnotVector degree=\" ")
        write(g,string(gs.degree))
        write(g,"\"> ")
        #    for i in 1:length(gs.knots)        write(g,string(gs.knots[i]," "))        end
        write(g,string(gs.knots," </KnotVector> \n"))
        write(g,"</Basis>\n")
        write(g,"</Basis>\n")
        write(g,"<coefs geoDim=\"3\">\n")
        for j in 1:2*lsz(gs)
            for i in 1:2*lsz(gs)
                write(g,string(gs.ctrpoints[1,idx(gs,gs.mesh.faces[f],i,j)]," ",gs.ctrpoints[2,idx(gs,gs.mesh.faces[f],i,j)]," ",gs.ctrpoints[3,idx(gs,gs.mesh.faces[f],i,j)],"\n"))
            end
        end
        write(g,"</coefs>\n")
        write(g,"</Geometry>\n")
    end

    write(g,string("<MultiPatch parDim=\"2\" id=\" ",F,"\" >\n"))
    write(g, string("<patches type=\"id_range\">",0," ",F-1,"</patches>\n"))
    write(g, "</MultiPatch>\n")
    #return 0
end

function SurfaceToGismo(gs,fl, F::Int64, mode = "w")
    g= open(string(fl,".xml"),mode)
    if mode =="w"
        write(g,"<?xml version=\"1.0\" encoding=\" UTF-8\" ?>\n")
        write(g,"<xml>\n")
    end
    AddSurfaceToGismoFile(gs,g,F)
    if mode == "w"
        write(g,"</xml>")
    end
    close(g)

    #return 0
end


#=function AddSurfaceToGismoFile(gs,g)

    for f in 1:nbf(gs.mesh)
        write(g,string("<Geometry type=\"TensorBSpline2\" id=\" ",f-1,"\" >\n"))
        write(g,"<Basis type=\"TensorBSplineBasis2\">\n")
        write(g,"<Basis type=\"BSplineBasis\" index=\" 0\"> \n")
        write(g,"<KnotVector degree=\" ")
        write(g,string(gs.degree))
        write(g,"\"> ")
        #for i in 1:length(gs.knots)            write(g,string(gs.knots[i]," "))        end
        write(g,string(gs.knots," </KnotVector> \n"))
        write(g,"</Basis>\n")
        write(g,"<Basis type=\"BSplineBasis\" index=\" 1\" >\n")
        write(g,"<KnotVector degree=\" ")
        write(g,string(gs.degree))
        write(g,"\"> ")
        #    for i in 1:length(gs.knots)        write(g,string(gs.knots[i]," "))        end
        write(g,string(gs.knots," </KnotVector> \n"))
        write(g,"</Basis>\n")
        write(g,"</Basis>\n")
        write(g,"<coefs geoDim=\"3\">\n")
        for j in 1:2*lsz(gs)
            for i in 1:2*lsz(gs)
                write(g,string(gs.ctrpoints[1,idx(gs,gs.mesh.faces[f],i,j)]," ",gs.ctrpoints[2,idx(gs,gs.mesh.faces[f],i,j)]," ",gs.ctrpoints[3,idx(gs,gs.mesh.faces[f],i,j)],"\n"))
            end
        end
        write(g,"</coefs>\n")
        write(g,"</Geometry>\n")
    end
    return 0
end=#
#SurfaceToGismo(bspline(gs),gs,"param")

# gs: GSpline
# basis: matrix with basis functions
# fl: filename to write to
# todo: use sparse matrix for basis: https://docs.julialang.org/en/v1/stdlib/SparseArrays/#man-csc
function BasisToGismo(gs, basis, fl, F::Int64, knt, mode="w")
    g= open(string(fl,".xml"),mode)

    if mode == "w"
        write(g,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        write(g,"<xml>\n")
    end

    local Fnew=F+1;

    for f in 1:F
        write(g,string("<Basis type=\"TensorBSplineBasis2\" id=\" ",Fnew,"\" >\n"))
        write(g,"<Basis type=\"BSplineBasis\" index=\" 0\"> \n")
        write(g,"<KnotVector degree=\" ")
        m, d = dim_deg(knt)
        write(g,string(d))
        write(g,"\"> ")

        write(g,string(knt," </KnotVector> \n"))
        write(g,"</Basis>\n")
        write(g,"<Basis type=\"BSplineBasis\" index=\" 1\" >\n")
        write(g,"<KnotVector degree=\" ")
        write(g,string(d))
        write(g,"\"> ")

        write(g,string(knt," </KnotVector> \n"))
        write(g,"</Basis>\n")
        write(g,"</Basis>\n")
        Fnew+=1;
    end


    write(g,string("<MultiBasis parDim=\"2\" id=\" ",Fnew,"\" >\n"))
    write(g, string("<patches type=\"id_range\">",F+1," ",Fnew-1,"</patches>\n"))
    write(g, "</MultiBasis>\n")




    write(g,string("<SparseMatrix rows=\" ", size(basis,1), "\" cols=\" ", size(basis,2),"\" id=\" ",Fnew+1,"\" >\n"))

    row=findnz(basis)[1];
    col=findnz(basis)[2];
    val=findnz(basis)[3];


    for i in 1:length(row)
        write(g,string(row[i]-1," ",col[i]-1," ",val[i],"\n"))
        #write(g,(@sprintf "%0.17g %0.17g %0.17g \n" row[i]-1 col[i]-1 val[i])); #*"\n")
    end

    write(g,"</SparseMatrix>\n")

    if mode == "w"
        write(g,"</xml>")
    end
    close(g)
    return 0
end


function MatrixtoGismo(matrix, fl, mode="w")

    g= open(string(fl,".xml"),mode)

    if mode == "w"
        write(g,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        write(g,"<xml>\n")
    end
    write(g,string("<Matrix rows=\" ", size(matrix,1), "\" cols=\" ", size(matrix,2),"\" id=\"0\">\n"))

    for i in 1:size(matrix,1)
        write(g,string(matrix[i,:],"\n"))
    end

    write(g,"</Matrix>\n")
    if mode == "w"
        write(g,"</xml>")
    end
    close(g)
    return 0
end

#BasisToGismo(gs,"basisGismo")

#
# Does both surface and basis in one file
#
"""
This function produces the `.xml` file containing both the spline geometry and associated biquintic basis functions (stored as sparse matrix) that can be loaded in G+Smo scripts to be used in IsoGeometric Analysis simulations or in point cloud fitting problems. 

It may also creates `.xml` files containing only the geometry or the set of basis functions.

## Example

    using GSplines
    m = offdata("triangle_planar.off")
    s1 = g1surface(m)
    basis = g1basis(m)
    ToGismo(s1,basis,s1.knots,"filename")

"""
function ToGismo(gs, basis, knt, fl)

    global F=nbf(gs.mesh);

    g=open(string(fl,".xml"),"w")

    write(g,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    write(g,"<xml>\n")

    close(g)

    SurfaceToGismo(gs,fl,F,"a")

    BasisToGismo(gs,basis,fl,F,knt,"a")

    g=open(string(fl,".xml"),"a")

    write(g,"</xml>")

    close(g)

    return 0
end

#not used
function GControlPointToGismoFormat(gs,fl)
        M=gs.basis[1]
        for i in 2:length(gs.basis)
                M=hcat(M,gs.basis[i])
        end
        g= open(string(fl,".xml"),"w")
        v=M\gs.ctrpoints'
        write(g,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        write(g,"<xml>\n")
        write(g,string("<Matrix rows=\" ",size(v,1) ,"\" cols=\" ",size(v,2) ,"\" id=\"0\">"))
        for i in 1:size(v,1)
                write(g,"\n")
                for j in 1:size(v,2)
                        write(g,string(v[i,j]," "))
                end
        end

        write(g,"\n </Matrix>\n")
        write(g,"</xml>\n")
        close(g)
end
#GControlPointToGismoFormat(gs,"Gcontroles")

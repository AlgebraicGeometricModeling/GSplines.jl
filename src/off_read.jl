export offread
#----------------------------------------------------------------------
"""
Read an off file and ouput a mesh.
```
offread("file.off")
```
"""
function offread(file::String)
    io = open(file)
    m = mesh(Float64)

    counts=0;
    nbv=0;
    nbf=0;
    for txt in eachline(io)
        l = split(txt, " ");
        if l[1] == "OFF"
            continue
        end
        if nbv==0
            nbv=parse(Int64, split(split(l[1],"//")[1],"/")[1]  );
            nbf=parse(Int64, split(split(l[2],"//")[1],"/")[1]  );
            continue
        end
        if counts<nbv
            counts+=1;
            pt = Float64[]
            s = 0
            for i in 1:length(l)
                if length(l[i])>0 && !startswith(l[i]," ") && s<3
                    push!(pt,parse(Float64,l[i]))
                    s+=1
                end
            end
            push_vertex!(m,pt)
        else
            f = Int64[]
            for i in 2:length(l)
                if length(l[i])>0 && !startswith(l[i]," ")
                    push!(f, 1+parse(Int64, split(split(l[i],"//")[1],"/")[1] ))
                end
            end
            push_face!(m,f)
        end
    end
    close(io)
    m
end

#----------------------------------------------------------------------
export offdata
"""
Read the off data in the 'file' of the data repository.
```
m = offdata("cube.off")
```

"""
function offdata(file::String)
    m = offread(joinpath(G1S[:pkgdir],"data/off",file))

end

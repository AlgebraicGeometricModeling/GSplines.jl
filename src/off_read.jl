export offread
#----------------------------------------------------------------------
"""
Read an off file and ouput a mesh.
```
offread("file.off")
```
"""
function offread(file::String)
    if !endswith(file,".off")
        file = file*".off"
    end
    
    io = open(file)
    m = mesh(Float64)

    ctrv = 0
    ctrf = 0
    ctre = 0
    nbv = 0
    nbf = 0
    nbe = 0
    for txt in eachline(io)
        l = split(txt, " ");
        if l[1] == "OFF"
            continue
        end
        if nbv==0
            nbv=parse(Int64, split(split(l[1],"//")[1],"/")[1]  );
            nbf=parse(Int64, split(split(l[2],"//")[1],"/")[1]  );
            nbe=parse(Int64, split(split(l[3],"//")[1],"/")[1]  );
            continue
        end
        if ctrv<nbv
            ctrv+=1;
            pt = Float64[]
            s = 0
            for i in 1:length(l)
                if length(l[i])>0 && !startswith(l[i]," ") && s<3
                    push!(pt,parse(Float64,l[i]))
                    s+=1
                end
            end
            push_vertex!(m,pt)
        elseif ctrf<nbf
            f = Int64[]
            for i in 2:length(l)
                if length(l[i])>0 && !startswith(l[i]," ")
                    push!(f, 1+parse(Int64, split(split(l[i],"//")[1],"/")[1] ))
                end
            end
            ctrf+=1
            push_face!(m,f)
        else
            if ctre<nbe
                e = parse.(Int64, l).+1
                ctre+=1
                push_edge!(m,e)
            end
        end
    end
    close(io)
    return m
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

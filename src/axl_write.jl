
function Axl.axlprint(io::IO, gs::GSpline, idt::Int64=0)
    S = bspline(gs)
    for i in 1:length(S)
        S[i][:color] =  Color(rand(0:255),rand(0:255),rand(0:255))
        Axl.axlprint(io,S[i],idt)
    end
end



export axldata
"""
Read the off data in the 'file' of the data repository.
```
m = axldata("y1m1.axl")
```

"""
function axldata(file::String)
    m = axlread(joinpath(G1S[:pkgdir],"data/axl",file))

end

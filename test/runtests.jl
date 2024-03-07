for f in readdir(pwd())
    if startswith(f,"test") && endswith(f,".jl")
        @info "Reading "*f
        try
            include(f)
        catch
            @error "Problem with "*f
            continue 
        end
    end
end

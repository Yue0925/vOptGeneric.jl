mutable struct SCP
    n :: Int64
    m :: Int64
    cover::Vector{Vector{Int64}}
    c :: Vector{Int64}
    name :: String
end



function SCP(n :: Int64, m :: Int64)
    return SCP(n, m, Vector{Vector{Int64}}(),
                    zeros(Int64,n), ""
                )
end



function readInstances(fname::String)::SCP
    f = open(fname)

    # initializing 
    vec = parse.(Int64, split(readline(f), " ")[2:3]) 
    n = vec[2] ; m = vec[1]
    inst = SCP(n, m)
    inst.name = split(split(fname, "/")[end], ".")[1]

    # reading cost 
    idx = 1
    while true
        vec = parse.(Int64, split(readline(f), " ")[2:end-1])
        inst.c[idx:idx+length(vec)-1] = vec[:]
        idx += length(vec)
        if idx > n break end 
    end

    # reading rows 
    for i in 1:m 
        # for each row 
        cols = parse(Int64, split(readline(f), " ")[2])
        push!(inst.cover, zeros(Int64, cols))
        idx = 1

        while true 
            vec = parse.(Int64, split(readline(f), " ")[2:end-1])
            inst.cover[i][idx:idx+length(vec)-1] = vec[:]
            idx += length(vec)
            if idx > cols break end 
        end
    end

    close(f)
    return inst
end


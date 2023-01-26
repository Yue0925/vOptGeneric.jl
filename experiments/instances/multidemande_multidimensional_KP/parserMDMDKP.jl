mutable struct MDMDKP 
    n :: Int64
    m :: Int64
    k :: Int64
    q :: Int64
    A_inf :: Matrix{Int64}
    A_sup :: Matrix{Int64}
    b_inf :: Vector{Int64}
    b_sup :: Vector{Int64}
    c :: Vector{Int64}
    name :: String
end

function MDMDKP(n :: Int64, m :: Int64, k :: Int64)
    return MDMDKP(n, m, k, 0, 
                    zeros(Int64,m, n), zeros(Int64,m, n), zeros(Int64,m), zeros(Int64,m), 
                    zeros(Int64,n), ""
                )
end

function readInstances(fname::String)::Vector{MDMDKP}
    f = open(fname)

    K = parse(Int64, readline(f)) ; Q = 6

    instances = Vector{MDMDKP}()

    # forech instacne in this file 
    for k in 1:K
        prefixName = split(split(fname, "/")[end], ".")[1] * "_" * string(k)

        vec = split(readline(f)[1:end-1], " ") 
        n = parse(Int64, vec[1]) ; m_infeq = parse(Int64, vec[end])

        inst = MDMDKP(n, m_infeq, k)

        # constraint ≤
        for i in 1:m_infeq
            a = split(readline(f)[1:end-1], " ") 
            inst.A_inf[i, :] = parse.(Int64, a)[:]
        end

        b = split(readline(f)[1:end-1], " ") 
        inst.b_inf[:] = parse.(Int64, b)[:]

        # constraint ≥
        for i in 1:m_infeq
            a = split(readline(f)[1:end-1], " ") 
            inst.A_sup[i, :] = parse.(Int64, a)[:]
        end

        b = split(readline(f)[1:end-1], " ") 
        inst.b_sup[:] = parse.(Int64, b)[:]

        # different objective cost 
        for q in 1:Q
            inst.q = q 
            instName = prefixName * "_$q" 
            c = split(readline(f)[1:end-1], " ") 
            inst.c[:] = parse.(Int64, c)[:] ; inst.name = instName
        end
        push!(instances, inst)
    end

    close(f)
    return instances
end


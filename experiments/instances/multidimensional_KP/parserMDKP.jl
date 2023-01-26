mutable struct MDKP 
    n :: Int64
    m :: Int64
    k :: Int64
    A :: Matrix{Int64}
    b :: Vector{Int64}
    c :: Vector{Int64}
    name :: String
end


function MDKP(n :: Int64, m :: Int64, k :: Int64)
    return MDKP(  n, m, k,
                    zeros(Int64,m, n), zeros(Int64,m), zeros(Int64,n), ""
                )
end

function readMDKP(fname::String)::Vector{MDKP}
    f = open(fname)

    K = parse(Int64, readline(f)) 
    instances = Vector{MDKP}()

    # foreach problem 
    for k in 1:K
        # read n & m, initialize instance 
        ss = split(readline(f), " ")[2:end] ; if ss[end] == "" pop!(ss) end 
        vec = parse.(Int64, ss )
        inst = MDKP(vec[1], vec[2], k ) ; inst.name = split(split(fname, "/")[end], ".")[1] * "_" * string(k)

        # read c
        s = 1
        while true
            ss = split(readline(f), " ")[2:end] ; if ss[end] == "" pop!(ss) end 
            vec = parse.(Int64, ss ) ; e = length(vec) + s - 1
            inst.c[s:e] = vec[:]
            s = e + 1
            if s > inst.n break end 
        end
        
        # ∀ j ∈ [1:m]
        for i in 1:inst.m 
            s = 1
            while true
                ss = split(readline(f), " ")[2:end] ; if ss[end] == "" pop!(ss) end 
                vec = parse.(Int64, ss ); e = length(vec) + s - 1
                inst.A[i, s:e] = vec[:] ; s = e + 1
                if s > inst.n break end 
            end
        end

        # read b 
        s = 1
        while true
            ss = split(readline(f), " ")[2:end] ;  if ss[end] == "" pop!(ss) end 
            vec = parse.(Int64, ss ) ; e = length(vec) + s - 1
            inst.b[s:e] = vec[:] ; s = e + 1
            if s > inst.m break end 
        end
        push!(instances, inst )

    end
    close(f)
    return instances
end

# readMDKP("./ORLib/mknap1.txt")
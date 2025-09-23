mutable struct HKP 
    n :: Int64
    A :: Vector{Int64}
    b :: Int64
    c :: Vector{Int64}
    name :: String
end



function HKP(n :: Int64)
    return HKP( n, zeros(Int64, n), 0, zeros(Int64,n), "" )
end


function readHKP(fname::String)::Vector{HKP}
    f = open(fname) ; K = 100
    instances = Vector{HKP}()

    # foreach instance in this file 
    for _ in 1:K 
        name = readline(f) ; m = 1
        n = parse(Int64, split(readline(f), " ")[2] )
        inst = HKP(n) ; inst.name = name 
        inst.b = parse(Int64, split(readline(f), " ")[2] )

        readline(f) ; readline(f) 
        for i = 1:inst.n 
            v = parse.(Int64, split(readline(f), ",")[2:3] )
            inst.c[i] = v[1] ; inst.A[i] = v[2]
        end
        push!(instances, inst )
        readline(f) ; readline(f)
    end

    close(f)
    return instances
end




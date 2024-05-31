mutable struct Forget
    n :: Int64
    m :: Vector{Int64} # 0 : >= ;  1 : <= ; 2 : == ; 
    c1 :: Vector{Int64}
    c2 :: Vector{Int64}
    A0 :: Matrix{Int64}
    b0 :: Vector{Int64}
    A1 :: Matrix{Int64}
    b1 :: Vector{Int64}
    A2 :: Matrix{Int64}
    b2 :: Vector{Int64}
    name :: String
end

# function Forget(n :: Int64, m :: Int64)
#     return Forget(n, m, zeros(Int64, n), zeros(Int64, n) , 
#                     Dict(0 => zeros(Int64, 0, 0) , 1 => zeros(Int64, 0, 0) , 2 => zeros(Int64, 0, 0) ), 
#                     Dict(0 => Vector{Int64}() , 1 => Vector{Int64}(), 2 => Vector{Int64}()), ""
#             )
# end


function readForget(fname :: String) :: Forget
    f = open(fname)
    vec = split(readline(f), " ")

    n = parse(Int64, vec[1]) ; m = parse(Int64, vec[2])

    # inst = Forget(n, m) ;
    A = zeros(Int64, m, n) ; c1 = zeros(Int64, n) ; c2 = zeros(Int64, n)

    readline(f) ; readline(f); readline(f)

    vec =  split(readline(f), " ")
    for i in 1:n 
        c1[i] = parse(Int64, vec[i])
    end

    vec =  split(readline(f), " ")
    for i in 1:n 
        c2[i] = parse(Int64, vec[i])
    end
    readline(f) ; readline(f)

    for i in 1:m
        vec =  split(readline(f), " ")
        for j in 1:n
            A[i, j] = parse(Int64, vec[j])
        end
    end

    readline(f) 

    lines = Dict{Int64, Vector{Int64}}(0 => [0], 1 => [0] , 2 => [0]) ; b = zeros(Int64, m)
    
    for i in 1:m
        vec =  split(readline(f), " ")
        push!(lines[parse(Int64, vec[1])], i)
        b[i] = parse(Int64, vec[2])
    end

    println(lines)

    println(b)
    
    l = length(lines[0]) - 1 ; A0 = zeros(Int64, l, n) ; b0 = zeros(Int64, l)
    if l >0
        for i in 
            
        end
    end

    close(f)
end

readForget("./instance/Forget20-UFLP_5_3_1-1000_1-100_spheredown_1_1.raw")

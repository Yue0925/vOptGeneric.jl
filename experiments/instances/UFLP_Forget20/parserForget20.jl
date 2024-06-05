mutable struct Forget
    n :: Int64
    m :: Dict{Int64, Vector{Int64}} # 0 : >= ;  1 : <= ; 2 : == ; 
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

    # println("A : ", size(A, 1), size(A, 2))

    # println("b : ", size(b))

    # println("lines : ", lines)
    
    l = length(lines[0]) - 1 ; A0 = zeros(Int64, l, n) ; b0 = zeros(Int64, l)
    if l >0
        for i in 2:l+1
            b0[i-1] = b[lines[0][i]]
            for j in 1:n
                A0[i-1, j] = A[lines[0][i], j]
            end
        end
    end

    l = length(lines[1]) - 1 ; A1 = zeros(Int64, l, n) ; b1 = zeros(Int64, l)
    if l >0
        for i in 2:l+1
            b1[i-1] = b[lines[1][i]]
            for j in 1:n
                A1[i-1, j] = A[lines[1][i], j]
            end
        end
    end

    l = length(lines[2]) - 1 ; A2 = zeros(Int64, l, n) ; b2 = zeros(Int64, l)
    if l >0
        for i in 2:l+1
            b2[i-1] = b[lines[2][i]]
            for j in 1:n
                A2[i-1, j] = A[lines[2][i], j]
            end
        end
    end

    close(f)

    return Forget(n, lines, c1, c2, A0, b0, A1, b1, A2, b2, split(split(fname, "/")[end], ".")[1])
end



function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N, X_E; total_time=nothing, infos=nothing)

    fout = open(outputName, "w")
    println(fout, "vars = $vars ; constr = $constr ")
  
    if method != :dicho && method != :epsilon
        println(fout, infos)
    else
      println(fout, "total_times_used = $total_time")
    end
    println(fout, "size_Y_N = ", length(Y_N))
    println(fout, "Y_N = ", Y_N)
    # println(fout)
    # println(fout, "size_X_E = ", length(X_E))

    close(fout)
  
    # displayGraphics(fname,Y_N, outputName)
end


# println(readForget("./instance/Forget20-UFLP_5_3_1-1000_1-100_spheredown_1_1.raw") )

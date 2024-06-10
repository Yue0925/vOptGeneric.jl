mutable struct ForgetKP
    n :: Int64
    c1 :: Vector{Int64}
    c2 :: Vector{Int64}
    A :: Vector{Int64}
    b :: Int64
    name :: String
end


function readForgetKP(fname :: String) :: ForgetKP
    f = open(fname)
    vec = split(readline(f), " ")

    n = parse(Int64, vec[1]) ; m = parse(Int64, vec[2])
    if m >= 2
        @error "multiple KP \t $fname"
    end

    A = zeros(Int64, n) ; c1 = zeros(Int64, n) ; c2 = zeros(Int64, n)

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

    # read ctr
    vec =  split(readline(f), " ")
    for j in 1:n
        A[j] = parse(Int64, vec[j])
    end

    readline(f) 

    # read b 
    vec =  split(readline(f), " ")
    if parse(Int64, vec[1]) != 1 @error "not <= ctr ! " end 

    b = parse(Int64, vec[2])

    close(f)

    return ForgetKP(n, c1, c2, A, b, split(split(fname, "/")[end], ".")[1])

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



# println(readForgetKP(ARGS[1]))
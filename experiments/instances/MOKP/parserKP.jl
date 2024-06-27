mutable struct MOKP
    n :: Int64
    c1 :: Vector{Int64}
    c2 :: Vector{Int64}
    A :: Vector{Int64}
    b :: Int64
    name :: String
end


function readMOKP(fname :: String) :: MOKP
    f = open(fname)

    p = parse(Int64, readline(f)) ; n = parse(Int64, readline(f))
    b = parse(Int64, readline(f))

    A = zeros(Int64, n) ; c1 = zeros(Int64, n) ; c2 = zeros(Int64, n)

    vec =split( split( split(readline(f), "[[")[end], "],")[1], ",")
    for i in 1:n 
        c1[i] = parse(Int64, vec[i])
    end

    vec =split( split( split(readline(f), "[")[end], "],")[1], ",")
    for i in 1:n 
        c2[i] = parse(Int64, vec[i])
    end

    for _ in 1:p-2
        readline(f) 
    end

    # read ctr
    vec = split( split( split(readline(f), "[")[end], "]")[1], ",")
    for j in 1:n
        A[j] = parse(Int64, vec[j])
    end


    close(f)

    return MOKP(n, c1, c2, A, b, split(split(fname, "/")[end], ".")[1])

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



# println(readMOKP(ARGS[1]))
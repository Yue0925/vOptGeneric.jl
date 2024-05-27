mutable struct AP
    n :: Int64
    C1 :: Matrix{Int64}
    C2 :: Matrix{Int64}
    name :: String
end

function AP(n :: Int64)
    return AP(n, zeros(Int64, n, n) , zeros(Int64, n, n), "")
end

function readAP(fname :: String)::AP
    f = open(fname)
    
    readline(f) # p , nb obj 
    n = parse(Int64, readline(f)) 
    ap = AP(n)

    # C1 
    for i in 1:n
        coeff =  split(split( split(readline(f), "[")[end], "]")[1], "," )
        for j in 1:n
            ap.C1[i, j] = parse(Int64, coeff[j] )
        end
    end

    readline(f) # , 

    # C2
    for i in 1:n
        coeff = split(split(split(readline(f), "[")[end], "]")[1], ",") 
        for j in 1:n
            ap.C2[i, j] = parse(Int64, coeff[j] )
        end
    end

    close(f)

    ap.name =  split(split(fname, "/")[end], ".")[1]

    return ap
    
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


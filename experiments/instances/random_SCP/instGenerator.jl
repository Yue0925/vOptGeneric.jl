

"""
Coefficients generator proposed by Pedersen et al, 2008.
"""
function generateC2(c1::Vector{Int64})::Vector{Int64}
    mini = minimum(c1) ; maxi = maximum(c1)
    c2 = zeros(Int64, length(c1))

    for i in 1:length(c1)
        if c1[i] < (maxi-mini)/2
            c2[i] = rand((maxi-c1[i]+mini):(maxi))
        elseif c1[i] >= (maxi-mini)/2
            c2[i] = rand((mini):(maxi-c1[i]+mini))
        else
            error("Coefficient error c1[$i] = $(c1[i]), maximum = $maxi, minimum = $mini !")
        end
    end
    return c2
end



function SCPgenerator(n::Int64, m::Int64, iden::Int64)
    println("n = $n , m = $m , iden = $iden ")
    # todo : change the density
    density_min = round(Int64, n/10) ; density_max = round(Int64, n*3/10)
    c1 = zeros(Int64 , n) 
    covers = Vector{Set{Int64}}(undef, m) ; for i = 1:m covers[i] = Set{Int64}() end 

    # every object is covered
    for elt in 1:n 
        set = rand(1:m) ; covers[set] = Set(elt)
    end

    # generate sets 
    for i = 1:m 
        len = rand(density_min:density_max)
        while length(covers[i]) < len 
            push!(covers[i], rand(1:n))
        end
    end

    # generate objectives 
    coeff_min = round(Int64, n/10)
    c1 = [rand(coeff_min:n) for _ in 1:n ]
    c2 = generateC2(c1)

    # write / stock instances 
    name = "SCP_" * string(n) * "_" * string(m) * "_" * string(iden)
    
    folder = "./instances/"
    if !isdir(folder)
        mkdir(folder)
    end

    fout = open(folder * name, "w")
    println(fout, "inst_name = \"$name\" ")
    println(fout, "n = $n ") ; println(fout, "m = $m ")
    println(fout, "c1 = $c1 ") ; println(fout, "c2 = $c2 ")
    println(fout, "Cover = $(collect.(covers)) ")
    close(fout)

end

for n in range(20, 100, step=10)
    for _ in 1:5 
        m = rand(round(Int64, n/10):round(Int64, 3*n/10))
        for iden in 1:3  
            SCPgenerator(n, m, iden)
        end
    end
end

n = 20
A = [582, 194, 679, 485, 396, 873, 594, 264, 462, 330, 582, 388, 291, 132, 660, 528, 970, 330, 582, 462] 
b=970 
c=[114, 38, 133, 95, 612, 171, 918, 408, 714, 510, 114, 76, 57, 204, 1020, 816, 190, 510, 114, 714]
sol = abs.(
    [-0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, 1.0, 0.10101010101009897, 0.7979797979798008, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0]
    )
    i =5, j = 2

    D = [16, 4, 20, 12, 18, 13, 14, 9]

# i = 5; j = 15

# 14, 5, 7
# 1, 1,  0.7
sol2 = abs.([-0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, 0.3566017316017326, -0.0, 0.054112554112553335, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0]
)

for i in [5, 8]
    for j in 1:n 
        if j !=i && A[i] + A[j] <= b

            @info "i =$i, j = $j "
                        
            O = sort(collect(1:n), by=x -> A[x], rev = true)
            println(" O = $O")
            I = []; D = []
            println("sol", sol)
            for o in O
                println("o = $o")
                if o == i || o == j || sol[o]>0.0
                    continue
                end

                if A[i]+A[o] > b
                    push!(I, o)
                else
                    push!(D, o)
                end

                println("D = $D \n I = $I ")
                
                xD = 0.0
                for d in D xD += sol[d] end 

                wD = 0
                for d in D wD += A[d] end 



                wI = 0
                for t in I wI += A[t] end 

                println("sol[i] - sol[j] - xD  = ", sol[i] - sol[j] - xD )
                println("sum(A) - A[i] - wD - wI = ", sum(A) - A[i] - wD - wI)
                if sol[i] - sol[j] - xD > 0.0 && sum(A) - A[i] - wD - wI <= b
                        println("D = $D")
                        @info "found"
                end

            end

        end
        
    end
end



# D = [3, 15, 1, 11, 19, 16, 4, 9, 12, 10, 18, 13, 8, 2, 7] 
# I = Any[17, 6] 
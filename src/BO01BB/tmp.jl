function test()
    colored = [true, true, false, true, false, true, true, false, true, true]

    println("size colored = $(length(colored)) ", colored)
    idx_l = colored[1] ? 0 : -1 ; idx_r = -1
    println("début idx_l = $idx_l  ->  idx_r = $idx_r")

    for i = 1:length(colored)
        if colored[i] continue end 

        if idx_l==-1 && (i+1)≤ length(colored) && colored[i+1]
            idx_l = i ; continue
        end

        if idx_l != -1 && idx_r==-1
            idx_r = i #; newSols = local_dichotomy(pb, LBS, idx_l, idx_r, round_results, verbose; args...)
            println("find interval idx_l = $idx_l  ->  idx_r = $idx_r")
            if (i+1)≤ length(colored) && colored[i+1]
                idx_l = i ; idx_r = -1
            else
                idx_l = -1; idx_r = -1
            end
        end
    end
    if idx_l != -1 && idx_r == -1
        idx_r = length(colored)+1
        println("find interval idx_l = $idx_l  ->  idx_r = $idx_r")
        # local_dichotomy(pb, LBS, idx_l, idx_r, round_results, verbose; args...)
    end
end


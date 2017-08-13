function sample_from_bounded_power_law(beta, xmin, xmax)
    a = xmin^-beta
    b = xmax^-beta - a
    return floor((a + b * rand())^(-1/beta))
end

function sample_chromosome(cell::Cell)
    ccount = size(cell.chromosomes)[1]
    chromosome_idx::UInt32 = floor(1+ccount*rand())
    while cell.chromosomes[chromosome_idx].size <= 1000000
        chromosome_idx = floor(1+ccount*rand())
    end
    return chromosome_idx
end

if haskey(options, "-s")
    function print_oldseg(idx, old_seg)
        print("chr$idx:")
        for seg in old_seg
            print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
        end
    end

    function print_newseg(new_seg)
        print(" ->")
        for seg in new_seg
            print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
        end
        println("")
    end

    function print_msg(msg)
        println(msg)
    end
else
    function print_oldseg(idx, old_seg) end
    function print_newseg(new_seg) end
    function print_msg(msg) end
end

function check_op_psum(ops)
    sum = 0
    for (lim, op) in ops
        sum += lim
    end
    if sum < 0.999999 || sum > 1.000001
        println("Chromosomal operation probabilities sum to $sum, not 1. Scaling accordingly.")
    end
    return sum
end

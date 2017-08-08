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

for file in filter(x -> endswith(x, ".jl"), readdir("operations"))
    include("operations/$file")
end

# 12% of focal amplifications and 26% of focal deletions are telomere-bound according to Zack/Beroukhim 2013 (doi:10.1038/ng.2760),
# whereas focal amplifications and deletions are relatively balanced (median of 11 amplifications and 12 deletions per sample).
# Frequency for whole genome duplications in CRC is stated as .43/sample. Following approximates these relative frequencies.

# Zack/Beroukhim also mentions arm-length or longer SCNAs (median 3 amplifications and 5 deletions per sample) and copy-neutral
# LOH events (median 1 per sample) (TODO: should these be implemented separately as well?)
# For now just estimate frequency of whole-chromosome events to amount to .57 per sample with 3 duplications for 5 deletions.

const ops = [(.43 * 1/24, whole_genome_duplication),            # 0.018
             (3/8 * .57 * 1/24, whole_chromosome_duplication),  # 0.009
             (5/8 * .57 * 1/24, whole_chromosome_loss),         # 0.015
             (.88 * 11/24, chromosome_tandem_duplication),      # 0.403
             (.06 * 11/24, chromosome_start_insertion),         # 0.028
             (.06 * 11/24, chromosome_end_insertion),           # 0.028
             (.74 * 12/24, chromosome_mid_deletion),            # 0.370
             (.13 * 12/24, chromosome_start_deletion),          # 0.065
             (.13 * 12/24, chromosome_end_deletion)]            # 0.065

const op_psum = check_op_psum(ops)

function chromosome_mutation(cell::Cell)
    r = rand() * op_psum
    for (p, op) in ops
        if(r < p)
            (c, k, data) = op(cell)
            ml = copy(cell.mutation_list)
            push!(ml, Mutation(op, data))
            return Cell(c, k, ml)
        end
        r -= p
    end
end

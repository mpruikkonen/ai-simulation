function already_duplicated(cell::Cell)
    for m in cell.mutation_list
        if m.op == whole_genome_duplication
            return true
        end
    end
    return false
end

function whole_genome_duplication(cell::Cell)
    if !already_duplicated(cell)
        new_chromosomes= Array{Chromosome, 1}(0)
        mcount = 0
        for chromosome in cell.chromosomes
            push!(new_chromosomes, chromosome)
            push!(new_chromosomes, chromosome)
            for driver in chromosome.drivers
                mcount += driver.gene.suppressor == false ? 1 : -1
            end
        end
        return new_chromosomes, cell.k+mcount, true
    else
        return cell.chromosomes, cell.k, false
    end
end

function whole_genome_duplication(segments, first_dup)
    if first_dup
        print_msg("genome duplicated")
        cidx = 1
        ccount = size(segments)[1]
        while ccount > 0
            splice!(segments, cidx+1:cidx, [copy(segments[cidx])])
            cidx += 2
            ccount -= 1
        end
    else
        print_msg("genome reduplication omitted")
    end
end

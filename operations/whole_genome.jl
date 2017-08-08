function whole_genome_duplication(cell::Cell)
    new_chromosomes= Array{Chromosome, 1}(0)
    mcount = 0
    for chromosome in cell.chromosomes
        push!(new_chromosomes, chromosome)
        push!(new_chromosomes, chromosome)
        for driver in chromosome.drivers
            mcount += driver.gene.suppressor == false ? 1 : -1
        end
    end
    return new_chromosomes, cell.k+mcount, 0
end

function whole_genome_duplication(segments, data)
    print_msg("genome duplicated")
    cidx = 1
    ccount = size(segments)[1]
    while ccount > 0
        splice!(segments, cidx+1:cidx, [copy(segments[cidx])])
        cidx += 2
        ccount -= 1
    end
end

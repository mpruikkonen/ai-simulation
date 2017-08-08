function whole_chromosome_duplication(cell::Cell)
    cidx = sample_chromosome(cell)
    chromosome = cell.chromosomes[cidx]
    mcount = 0
    for driver in chromosome.drivers
        mcount += driver.gene.suppressor == false ? 1 : -1
    end
    new_chromosomes = copy(cell.chromosomes)
    splice!(new_chromosomes, cidx+1:cidx, [new_chromosomes[cidx]])
    return new_chromosomes, cell.k+mcount, cidx
end

function whole_chromosome_loss(cell::Cell)
    cidx = sample_chromosome(cell)
    chromosome = cell.chromosomes[cidx]
    mcount = 0
    for driver in chromosome.drivers
        mcount += driver.gene.suppressor == true ? 1 : -1
    end
    new_chromosomes = copy(cell.chromosomes)
    splice!(new_chromosomes, cidx)
    return new_chromosomes, cell.k+mcount, cidx
end

function whole_chromosome_loss(segments, cidx)
    print_msg("chr$cidx deleted")
    splice!(segments, cidx)
end

function whole_chromosome_duplication(segments, cidx)
    print_msg("chr$cidx duplicated")
    splice!(segments, cidx+1:cidx, [copy(segments[cidx])])
end

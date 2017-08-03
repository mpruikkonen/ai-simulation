# According to the description for the supplementary figure 2a in Zack/Beroukhim 2013 (doi:10.1038/ng.2760)
# chromosomal event frequencies for amplifications of more than 400kb and less than a chromosome arm in length
# follow a power law f(L) = 1 / L^b where b=0.45 for amplifications reaching the telomere and b=1.12 for
# internal events. This file implements chromosomal amplifications accordingly.

function chromosome_start_insertion(cell::Cell)
    cidx = sample_chromosome(cell)
    chromosome = cell.chromosomes[cidx]
    ccount = size(cell.chromosomes)[1]
    fidx::UInt32 = floor(1 + ccount * rand())
    while cell.chromosomes[fidx].centromere_pos <= 400000
        fidx = floor(1 + ccount * rand())
    end
    from_chromosome = cell.chromosomes[fidx]
    # we only sample from beginnings of chromosomes for now for start_insertions and from ends for end_insertions
    # TODO insertions with inverted sequences?
    insert_bases::UInt32 = sample_from_bounded_power_law(0.45, 400000, from_chromosome.centromere_pos-1)

    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        push!(new_drivers, Driver_Locus(driver.pos+insert_bases, driver.gene))
    end
    for driver in from_chromosome.drivers
        if driver.pos < insert_bases
            push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
            mcount += driver.gene.suppressor == false ? 1 : -1
        end
    end
    new_chromosomes = copy(cell.chromosomes)
    new_chromosomes[cidx] = Chromosome(chromosome.size+insert_bases, chromosome.centromere_pos+insert_bases, new_drivers)
    return new_chromosomes, cell.k+mcount, (cidx, fidx, insert_bases)
end

function chromosome_end_insertion(cell::Cell)
    cidx = sample_chromosome(cell)
    chromosome = cell.chromosomes[cidx]
    ccount = size(cell.chromosomes)[1]
    fidx::UInt32 = floor(1 + ccount * rand())
    while cell.chromosomes[fidx].size - cell.chromosomes[fidx].centromere_pos <= 400000
        fidx = floor(1 + ccount * rand())
    end
    from_chromosome = cell.chromosomes[fidx]
    insert_bases::UInt32 = sample_from_bounded_power_law(0.45, 400000, from_chromosome.size-from_chromosome.centromere_pos)

    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
    end
    for driver in from_chromosome.drivers
        if driver.pos < insert_bases
            push!(new_drivers, Driver_Locus(driver.pos-(from_chromosome.size-insert_bases)+chromosome.size, driver.gene))
            mcount += driver.gene.suppressor == false ? 1 : -1
        end
    end
    new_chromosomes = copy(cell.chromosomes)
    new_chromosomes[cidx] = Chromosome(chromosome.size+insert_bases, chromosome.centromere_pos, new_drivers)
    return new_chromosomes, cell.k+mcount, (cidx, fidx, insert_bases)
end

function chromosome_tandem_duplication(cell::Cell)
    cidx = sample_chromosome(cell)
    chromosome = cell.chromosomes[cidx]
    break_len::UInt32 = break1::UInt32 = break2::UInt32 = 0
    if chromosome.centromere_pos > 500000 && rand() < chromosome.centromere_pos / chromosome.size
        break_len = sample_from_bounded_power_law(1.12, 400000, chromosome.centromere_pos-1)
        break1 = floor(1 + rand() * (chromosome.centromere_pos - 1 - break_len))
        break2 = break1 + break_len
    else
        break_len = sample_from_bounded_power_law(1.12, 400000, chromosome.size-chromosome.centromere_pos)
        break1 = floor(chromosome.centromere_pos + rand() * ((chromosome.size - chromosome.centromere_pos) - break_len))
        break2 = break1 + break_len
    end

    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
        if(driver.pos > break1 && driver.pos < break2)
            push!(new_drivers, Driver_Locus(driver.pos+break_len, driver.gene))
            mcount += driver.gene.suppressor == false ? 1 : -1
        end
    end
    new_chromosomes = copy(cell.chromosomes)
    if break1 < chromosome.centromere_pos
        new_chromosomes[cidx] = Chromosome(chromosome.size+break_len, chromosome.centromere_pos+break_len, new_drivers)
    else
        new_chromosomes[cidx] = Chromosome(chromosome.size+break_len, chromosome.centromere_pos, new_drivers)
    end
    return new_chromosomes, cell.k+mcount, (cidx, break1, break2)
end

##############################################################################
# Reference sequence segment generation functions for reporting / statistics #
##############################################################################
function chromosome_start_insertion(segments, data)
    (cidx, fidx, bp) = data
    print_oldseg(cidx, segments[cidx])
    i = 1
    while bp > 0
        if bp >= segments[fidx][i][3] - segments[fidx][i][2]
            bp -= segments[fidx][i][3] - segments[fidx][i][2]
            splice!(segments[cidx], i:i-1, [segments[fidx][i]])
            i += 1
        else
            splice!(segments[cidx], i:i-1, [(segments[fidx][i][1], segments[fidx][i][2], segments[fidx][i][2]+bp)])
            bp = 0
        end
    end
    print_newseg(segments[cidx])
end

function chromosome_end_insertion(segments, data)
    (cidx, fidx, bp) = data
    print_oldseg(cidx, segments[cidx])
    fi = size(segments[fidx])[1]
    ci = size(segments[cidx])[1]
    while bp > 0
        if bp >= segments[fidx][fi][3] - segments[fidx][fi][2]
            bp -= segments[fidx][fi][3] - segments[fidx][fi][2]
            splice!(segments[cidx], ci+1:ci, [segments[fidx][fi]])
            fi -= 1
        else
            splice!(segments[cidx], ci+1:ci, [(segments[fidx][fi][1], segments[fidx][fi][3]-bp, segments[fidx][fi][3])])
            bp = 0
        end
    end
    print_newseg(segments[cidx])
end

function duplicate_inserted_segments(new_segments, segments, bp, bp2)
    i = 1
    while true
        if bp >= segments[i][3] - segments[i][2]
            bp -= segments[i][3] - segments[i][2]
            bp2 -= segments[i][3] - segments[i][2]
        else
            while bp2 >= segments[i][3] - segments[i][2]
                bp2 -= segments[i][3] - segments[i][2]
                push!(new_segments, (segments[i][1], segments[i][2]+bp, segments[i][3]))
                bp = 0
                i += 1
            end
            if bp2 > 0
                push!(new_segments, (segments[i][1], segments[i][2]+bp, segments[i][2]+bp2))
            end
            return
        end
        i += 1
    end
end

function chromosome_tandem_duplication(chr_segments, data)
    cidx = data[1]
    segments = chr_segments[cidx]
    print_oldseg(cidx, segments)
    new_segments = Array{Tuple{String, UInt32, UInt32}, 1}()
    bp2 = data[3]
    i = 1
    while true
        if bp2 >= segments[i][3] - segments[i][2]
            push!(new_segments, segments[i])
            bp2 -= segments[i][3] - segments[i][2]
        else
            if bp2 > 0
                push!(new_segments, (segments[i][1], segments[i][2], segments[i][2]+bp2))
            end
            duplicate_inserted_segments(new_segments, segments, data[2], data[3])
            push!(new_segments, (segments[i][1], segments[i][2]+bp2, segments[i][3]))
            i += 1
            while(i <= size(segments)[1])
                push!(new_segments, segments[i])
                i += 1
            end
            break
        end
        i += 1
    end
    chr_segments[cidx] = new_segments
    print_newseg(new_segments)
end

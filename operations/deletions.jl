# According to the description for the supplementary figure 2b in Zack/Beroukhim 2013 (doi:10.1038/ng.2760)
# chromosomal event frequencies for deletions of more than 400kb and less than a chromosome arm in length
# follow a power law f(L) = 1 / L^b where b=0.36 for end deletions reaching the telomere and b=1.15 for
# internal (double-break) events. This file implements chromosomal deletions accordingly.

function sample_from_bounded_power_law(beta, xmin, xmax)
    a = xmin^-beta
    b = xmax^-beta - a
    return floor((a + b * rand())^(-1/beta))
end

function get_deletion_mcount(cell::Cell, gene::Driver_Gene)
    gcount = 0
    for chromosome in cell.chromosomes
        for driver in chromosome.drivers
            if driver.gene == gene
                gcount += 1
            end
        end
    end

    if gene.suppressor == true
        return gcount == 1 ? 1 : 0 # fitness increases when losing last copy of suppressor
    else
        return gcount > 2 ? -1 : 0 # fitness decreases when losing extra copies of oncogene
    end
end

function chromosome_start_deletion(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    delete_bases::UInt32 = sample_from_bounded_power_law(0.36, 400000, chromosome.centromere_pos)
    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        if(driver.pos > delete_bases)
            push!(new_drivers, Driver_Locus(driver.pos-delete_bases, driver.gene))
        else
            mcount += get_deletion_mcount(cell, driver.gene)
        end
    end
    push_c1!(Chromosome(chromosome.size-delete_bases, chromosome.centromere_pos-delete_bases, new_drivers))
    push_c2!(chromosome)
    return mcount, 0, delete_bases
end

function chromosome_end_deletion(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    delete_bases::UInt32 = sample_from_bounded_power_law(0.36, 400000, chromosome.size-chromosome.centromere_pos)
    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        if(driver.pos < chromosome.size-delete_bases)
            push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
        else
            mcount += get_deletion_mcount(cell, driver.gene)
        end
    end
    push_c1!(Chromosome(chromosome.size-delete_bases, chromosome.centromere_pos, new_drivers))
    push_c2!(chromosome)
    return mcount, 0, delete_bases
end

function chromosome_mid_deletion(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    break_len::UInt32 = break1::UInt32 = break2::UInt32 = 0
    if rand() < chromosome.centromere_pos / chromosome.size
        break_len = sample_from_bounded_power_law(1.15, 400000, chromosome.centromere_pos)
        break1 = floor(rand() * (chromosome.centromere_pos - break_len))
        break2 = break1 + break_len
    else
        break_len = sample_from_bounded_power_law(1.15, 400000, chromosome.size-chromosome.centromere_pos)
        break1 = floor(chromosome.centromere_pos + rand() * ((chromosome.size - chromosome.centromere_pos) - break_len))
        break2 = break1 + break_len
    end
            
    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        if(driver.pos < break1 || driver.pos > break2)
            push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
        else
            mcount += get_deletion_mcount(cell, driver.gene)
        end
    end
    if break1 < chromosome.centromere_pos
        push_c1!(Chromosome(chromosome.size-break_len, chromosome.centromere_pos-break_len, new_drivers))
    else
        push_c1!(Chromosome(chromosome.size-break_len, chromosome.centromere_pos, new_drivers))
    end
    push_c2!(chromosome)
    return mcount, 0, (break1, break2)
end

##############################################################################
# Reference sequence segment generation functions for reporting / statistics #
##############################################################################
function chromosome_start_deletion(segments::Array{Tuple{String, UInt32, UInt32}, 1}, cell1::Bool, bp::UInt32, push_seg!::Function)
    if cell1 == true
        while bp > 0
            if bp >= segments[1][3] - segments[1][2]
                bp -= segments[1][3] - segments[1][2]
                shift!(segments)
            else
                segments[1] = (segments[1][1], segments[1][2]+bp, segments[1][3])
                bp = 0
            end
        end
    end
    push_seg!(segments)
end

function chromosome_end_deletion(segments::Array{Tuple{String, UInt32, UInt32}, 1}, cell1::Bool, bp::UInt32, push_seg!::Function)
    if cell1 == true
        while bp > 0
            last = size(segments)[1]
            if bp >= segments[last][3] - segments[last][2]
                bp -= segments[last][3] - segments[last][2]
                pop!(segments)
            else
                segments[last] = (segments[last][1], segments[last][2], segments[last][3]-bp)
                bp = 0
            end
        end
    end
    push_seg!(segments)
end

function chromosome_mid_deletion(segments::Array{Tuple{String, UInt32, UInt32}, 1}, cell1::Bool, breaks::Tuple{UInt32, UInt32}, push_seg!::Function)
    if cell1 == true
        new_segments = Array{Tuple{String, UInt32, UInt32}, 1}()
        bp = breaks[1]
        i = 1
        while true
            if bp >= segments[i][3] - segments[i][2]
                push!(new_segments, segments[i])
                bp -= segments[i][3] - segments[i][2]
            else
                push!(new_segments, (segments[i][1], segments[i][2], segments[i][2]+bp))
                deleted_bases = breaks[2]-breaks[1]
                if(deleted_bases < (segments[i][3] - (segments[i][2]+bp)))
                    push!(new_segments, (segments[i][1], segments[i][2]+bp+deleted_bases, segments[i][3]))
                else
                    deleted_bases -= segments[i][3] - (segments[i][2]+bp)
                    i += 1
                    while deleted_bases >= segments[i][3] - segments[i][2]
                        deleted_bases -= segments[i][3] - segments[i][2]
                        i += 1
                    end
                    push!(new_segments, (segments[i][1], segments[i][2]+deleted_bases, segments[i][3]))
                end
                i += 1
                while(i <= size(segments)[1])
                    push!(new_segments, segments[i])
                    i += 1
                end
                break
            end
            i += 1
        end
        push_seg!(new_segments)
    else
        push_seg!(segments)
    end
end

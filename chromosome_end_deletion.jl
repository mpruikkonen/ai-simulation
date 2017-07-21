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
    # TODO use data of actual centromere position?
    delete_bases::UInt32 = floor(rand()*(chromosome.size/2))
    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        if(driver.pos > delete_bases)
            push!(new_drivers, Driver_Locus(driver.pos-delete_bases, driver.gene))
        else
            mcount += get_deletion_mcount(cell, driver.gene)
        end
    end
    push_c1!(Chromosome(chromosome.size-delete_bases, new_drivers))
    push_c2!(chromosome)
    return mcount, 0, delete_bases
end

function chromosome_end_deletion(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    # TODO use data of actual centromere position?
    delete_bases::UInt32 = floor(rand()*(chromosome.size/2))
    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    for driver in chromosome.drivers
        if(driver.pos < delete_bases)
            push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
        else
            mcount += get_deletion_mcount(cell, driver.gene)
        end
    end
    push_c1!(Chromosome(chromosome.size-delete_bases, new_drivers))
    push_c2!(chromosome)
    return mcount, 0, delete_bases
end

function chromosome_start_deletion(segments::Array{Tuple{String, UInt32, UInt32}, 1}, cell1::Bool, bp::UInt32, push_seg!::Function)
    if cell1 == true
        while bp > 0
            if bp >= segments[1][3] - segments[1][2]
                print("start: $bp >= $(segments[1][3]) - $(segments[1][2])")
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
                print("end: $bp >= $(segments[1][3]) - $(segments[1][2])")
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


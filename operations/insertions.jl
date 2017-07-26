# According to the description for the supplementary figure 2a in Zack/Beroukhim 2013 (doi:10.1038/ng.2760)
# chromosomal event frequencies for amplifications of more than 400kb and less than a chromosome arm in length
# follow a power law f(L) = 1 / L^b where b=0.45 for amplifications reaching the telomere and b=1.12 for
# internal events. This file implements chromosomal amplifications accordingly.

function sample_from_bounded_power_law(beta, xmin, xmax)
    a = xmin^-beta
    b = xmax^-beta - a
    return floor((a + b * rand())^(-1/beta))
end

function chromosome_tandem_duplication(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    break_len::UInt32 = break1::UInt32 = break2::UInt32 = 0
    if rand() < chromosome.centromere_pos / chromosome.size
        break_len = sample_from_bounded_power_law(1.12, 400000, chromosome.centromere_pos)
        break1 = floor(rand() * (chromosome.centromere_pos - break_len))
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
            if driver.gene.suppressor == false
               mcount += 1 # TODO do we want just the first duplication of an oncogene to increase fitness instead of each one?
            end
        end
    end
    if break1 < chromosome.centromere_pos
        push_c1!(Chromosome(chromosome.size+break_len, chromosome.centromere_pos+break_len, new_drivers))
    else
        push_c1!(Chromosome(chromosome.size+break_len, chromosome.centromere_pos, new_drivers))
    end
    push_c2!(chromosome)
    return mcount, 0, (break1, break2)
end

##############################################################################
# Reference sequence segment generation functions for reporting / statistics #
##############################################################################
function duplicate_inserted_segments(new_segments::Array{Tuple{String, UInt32, UInt32}, 1}, segments::Array{Tuple{String, UInt32, UInt32}, 1}, bp::UInt32, bp2::UInt32)
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

function chromosome_tandem_duplication(segments::Array{Tuple{String, UInt32, UInt32}, 1}, cell1::Bool, breaks::Tuple{UInt32, UInt32}, push_seg!::Function)
    if cell1 == true
        new_segments = Array{Tuple{String, UInt32, UInt32}, 1}()
        bp2 = breaks[2]
        i = 1
        while true
            if bp2 >= segments[i][3] - segments[i][2]
                push!(new_segments, segments[i])
                bp2 -= segments[i][3] - segments[i][2]
            else
                if bp2 > 0
                    push!(new_segments, (segments[i][1], segments[i][2], segments[i][2]+bp2))
                end
                duplicate_inserted_segments(new_segments, segments, breaks[1], breaks[2])
                if segments[i][3] > segments[i][2]+bp2
                    push!(new_segments, (segments[i][1], segments[i][2]+bp2, segments[i][3]))
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

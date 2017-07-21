function whole_chromosome_duplication(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    push_c1!(chromosome)
    push_c1!(chromosome)
    push_c2!(chromosome)
    mcount = 0
    for driver in chromosome.drivers
        if driver.gene.suppressor == false
            mcount += 1 # TODO do we want to count just the first duplication of an oncogene instead of each one?
        end
    end
    return mcount, 0, 0
end

function whole_chromosome_duplication(segments::Array{Tuple{String, UInt32, UInt32}, 1}, cell1::Bool, data, push_seg!::Function)
    push_seg!(segments)
    if cell1 == true
        push_seg!(segments)
    end
end


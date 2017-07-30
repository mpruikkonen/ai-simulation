function write_bed(regions)
    fname = options["-o"]
    open(fname, "w") do f
        if haskey(options, "--std-bed")
            for (key, val) in regions
                write(f, "chr$(key[1])\t$(key[2])\t$(key[3])\t\t$(val[1])\n")
            end
        else
            for (key, val) in regions
                write(f, "chr$(key[1])\t$(key[2])\t$(key[3])\t\t$(val[1])\t$(val[2])\t$(val[3])\n")
            end
        end
    end
    println("Wrote BED file to '$fname'")
end

function get_region_name(chr, position)
    for (pos, loc) in hg38_cytobands[chr]
        if(pos >= position)
            return loc
        end
    end
end

function print_chromosome_regions(segments)
    regions = Dict{Tuple{String, Int, Int}, Tuple{Int, Int, Int}}()
    println("--")
    println("Final karyotype of the cell with the highest fitness:")
    idx = 1
    for chr in segments
        print("chr$idx:")
        scount = 0
        for seg in chr
            print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
            scount += 1
        end
        print(" (")
        for seg in chr
            c = seg[1][2:end]
            tcount = fcount = mcount = 0
            if haskey(regions, ("$c", seg[2], seg[3]))
                (tcount, fcount, mcount) = regions[("$c", seg[2], seg[3])]
            end
            tcount += 1
            if(seg[1][1] == 'f')
                fcount += 1
            else
                mcount += 1
            end
            regions[("$c", seg[2], seg[3])] = (tcount, fcount, mcount)

            b = get_region_name(c, seg[2])
            e = get_region_name(c, seg[3])
            scount -= 1
            if scount == 0
                println("$c$b->$c$e)")
            else
                print("$c$b->$c$e ")
            end
        end
        idx += 1
    end
    write_bed(regions)
end

function print_mutation(idx, old_seg, new_seg)
    print("chr$idx:")
    for seg in old_seg
        print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
    end
    print(" ->")
    for seg in new_seg
        print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
    end
end

function get_segments_for_clean_cell()
    segments = Array{Array{Tuple{String, UInt32, UInt32}, 1}, 1}()
    for i in 1:22
        push!(segments, [("m$i", 0, hg38_chromosome_sizes["$i"])])
        push!(segments, [("f$i", 0, hg38_chromosome_sizes["$i"])])
    end
    push!(segments, [("mX", 0, hg38_chromosome_sizes["X"])])
    push!(segments, [("fY", 0, hg38_chromosome_sizes["Y"])])
    return segments
end

function print_cell(cell::Cell)
    segments = get_segments_for_clean_cell()

    for m in cell.mutation_list
        push_count = 0
        printed_something = false
        function push_seg!(s::Array{Tuple{String, UInt32, UInt32}, 1})
            if push_count == 0
                if s != segments[m.chromosome_idx]
                    print_mutation(m.chromosome_idx, segments[m.chromosome_idx], s)
                    printed_something = true
                end
                segments[m.chromosome_idx] = s
            else
                if printed_something == false
                    print_mutation(m.chromosome_idx, segments[m.chromosome_idx], segments[m.chromosome_idx])
                    printed_something = true
                end
                print(",")
                for seg in s
                    print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
                end
                splice!(segments, m.chromosome_idx+1:m.chromosome_idx, [copy(s)])
            end
            push_count += 1
        end

        # Dispatch to chromosomal operation -specific refseq segment generator function, which should call push_seq! with
        # an array of ("reference chromosome name", start_bp_of_segment, end_bp_of_segment) tuples for each chromosome
        # that was created. Second parameter indicates whether this was the first or second output cell of the chromosomal
        # mutation operation, and the third parameter contains any data that was returned by the mutation operation.
        m.op(copy(segments[m.chromosome_idx]), m.cell1, m.data, push_seg!)

        if push_count == 0
            print("Chromosome $(m.chromosome_idx) deleted")
            printed_something = true
            splice!(segments, m.chromosome_idx:m.chromosome_idx)
        end
        if printed_something == true
            println("")
        end
    end
    
    print_chromosome_regions(segments)
end


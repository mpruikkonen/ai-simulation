function write_bed_file(chr_regions)
    fname = options["-o"]
    open(fname, "w") do f
        if haskey(options, "--std-bed")
            for chr in sort(collect(keys(chr_regions)))
                for r in chr_regions[chr]
                    write(f, "chr$chr\t$(r[1])\t$(r[2])\t\t$(r[3])\n")
                end
            end
        else
            for chr in sort(collect(keys(chr_regions)))
                for r in chr_regions[chr]
                    write(f, "chr$chr\t$(r[1])\t$(r[2])\t\t$(r[3])\t$(r[4])\t$(r[5])\n")
                end
            end
        end
    end
    println("Wrote BED file to '$fname'")
end

function add_seg_end!(regions, idx, cur_pos, seg_end, p)
    while idx <= size(regions)[1]
        if regions[idx][1] > cur_pos
            if regions[idx][1] >= seg_end
                splice!(regions, idx:idx-1, [(cur_pos, seg_end, 1, p == 'f' ? 1 : 0, p == 'm' ? 1 : 0)])
                return
            else
                splice!(regions, idx:idx-1, [(cur_pos, regions[idx][1], 1, p == 'f' ? 1 : 0, p == 'm' ? 1 : 0)])
                idx += 1
                cur_pos = regions[idx][1]
            end
        end
        if regions[idx][2] >= seg_end
            if regions[idx][2] > seg_end
                splice!(regions, idx+1:idx, [(seg_end, regions[idx][2], regions[idx][3], regions[idx][4], regions[idx][5])])
            end
            regions[idx] = (cur_pos, seg_end, regions[idx][3]+1, regions[idx][4] + (p == 'f' ? 1 : 0), regions[idx][5] + (p == 'm' ? 1 : 0))
            return
        else
            regions[idx] = (cur_pos, regions[idx][2], regions[idx][3]+1, regions[idx][4] + (p == 'f' ? 1 : 0), regions[idx][5] + (p == 'm' ? 1 : 0))
            cur_pos = regions[idx][2]
        end
        idx += 1
    end
    push!(regions, (cur_pos, seg_end, 1, p == 'f' ? 1 : 0, p == 'm' ? 1 : 0))
end

function add_seg!(regions, seg_start, seg_end, p)
    idx = 1
    while idx <= size(regions)[1]
        if regions[idx][2] > seg_start
            if regions[idx][1] <= seg_start
                if regions[idx][1] < seg_start
                    splice!(regions, idx:idx-1, [(regions[idx][1], seg_start, regions[idx][3], regions[idx][4], regions[idx][5])])
                    idx += 1
                end
                if regions[idx][2] > seg_end
                    splice!(regions, idx+1:idx, [(seg_end, regions[idx][2], regions[idx][3], regions[idx][4], regions[idx][5])])
                    regions[idx] = (seg_start, seg_end, regions[idx][3]+1, regions[idx][4] + (p == 'f' ? 1 : 0), regions[idx][5] + (p == 'm' ? 1 : 0))
                else
                    regions[idx] = (seg_start, regions[idx][2], regions[idx][3]+1, regions[idx][4] + (p == 'f' ? 1 : 0), regions[idx][5] + (p == 'm' ? 1 : 0))
                    if regions[idx][2] < seg_end
                        add_seg_end!(regions, idx+1, regions[idx][2], seg_end, p)
                    end
                end
            else
                if regions[idx][1] >= seg_end
                    splice!(regions, idx:idx-1, [(seg_start, seg_end, 1, p == 'f' ? 1 : 0, p == 'm' ? 1 : 0)])
                else
                    splice!(regions, idx:idx-1, [(seg_start, regions[idx][1], 1, p == 'f' ? 1 : 0, p == 'm' ? 1 : 0)])
                    add_seg_end!(regions, idx+1, regions[idx+1][1], seg_end, p)
                end
            end
            return
        end
        idx += 1
    end
    push!(regions, (seg_start, seg_end, 1, p == 'f' ? 1 : 0, p == 'm' ? 1 : 0))
end

function get_region_name(chr, position)
    for (pos, loc) in hg38_cytobands[chr]
        if(pos >= position)
            return loc
        end
    end
end

function print_chromosome_regions(segments)
    chr_regions = Dict{String, Array{Tuple{Int, Int, Int, Int, Int}}}()
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
            parent = seg[1][1]
            if haskey(chr_regions, "$c")
                add_seg!(chr_regions["$c"], seg[2], seg[3], parent)
            else
                chr_regions["$c"] = [(seg[2], seg[3], 1, parent == 'f' ? 1 : 0, parent == 'm' ? 1 : 0)]
            end
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
    write_bed_file(chr_regions)
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


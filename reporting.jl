function write_bedgraph_file(chr_regions)
    fname = options["-o"]
    open(fname, "w") do f
        if haskey(options, "--std-bed")
            for chr in sort(collect(keys(chr_regions)))
                for r in chr_regions[chr]
                    write(f, "$chr\t$(r[1])\t$(r[2])\t$(r[3])\n")
                end
            end
        else
            for chr in sort(collect(keys(chr_regions)))
                for r in chr_regions[chr]
                    write(f, "$chr\t$(r[1])\t$(r[2])\t$(r[3])\t$(r[4])\t$(r[5])\n")
                end
            end
        end
    end
    if haskey(options, "-b")
        println("""Wrote BedGraph file to '$fname' and Bed files for each iteration to '$(options["-b"])'""")
    else
        println("Wrote BedGraph file to '$fname'")
    end
end

function write_segments_to_bed_file(segments, iteration)
    fname = """$(options["-b"])/$(iteration).bed"""
    open(fname, "w") do f
        for chr in segments
            for seg in chr
                c = seg[1][2:end]
                write(f, "$c\t$(seg[2])\t$(seg[3])\n")
            end
        end
    end
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
    for (pos, loc) in ref_cytobands[chr]
        if(pos >= position)
            return loc
        end
    end
end

function add_segments_to_chr_regions(segments, chr_regions)
    idx = 1
    for chr in segments
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
        end
        idx += 1
    end
end

function get_segments_for_clean_cell()
    segments = Array{Array{Tuple{String, UInt32, UInt32}, 1}, 1}()
    for i in 1:22
        push!(segments, [("m$i", 0, ref_chromosome_sizes["$i"])])
        push!(segments, [("f$i", 0, ref_chromosome_sizes["$i"])])
    end
    push!(segments, [("mX", 0, ref_chromosome_sizes["X"])])
    push!(segments, [("fY", 0, ref_chromosome_sizes["Y"])])
    return segments
end

function print_cell(cell::Cell, chr_regions, iteration)
    segments = get_segments_for_clean_cell()
    for m in cell.mutation_list
        # Dispatch to chromosomal operation -specific refseq segment generator function, which should alter the provided segment
        # array of ("reference chromosome name", start_bp_of_segment, end_bp_of_segment) tuples according to the mutation operation
        # that was performed. Second parameter contains any data that was returned by the mutation operation.
        m.op(segments, m.data)
    end
    add_segments_to_chr_regions(segments, chr_regions)
    if haskey(options, "-b")
        write_segments_to_bed_file(segments, iteration)
    end
end


function write_bedgraph_file(fname, chr_regions)
    open(fname, "w") do f
        for chr in sort(collect(keys(chr_regions)))
            for r in chr_regions[chr]
                write(f, "$chr\t$(r[1])\t$(r[2])\t$(r[3])\n")
            end
        end
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

function add_segments_to_region_array(segments, chr_regions)
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
        end
        idx += 1
    end
end

function add_segments_to_ai_arrays(segments, ai_gains, ai_losses)
    chr_regions = Dict{String, Array{Tuple{Int, Int, Int, Int, Int}}}()
    add_segments_to_region_array(segments, chr_regions)
    for c in sort(collect(keys(chr_regions)))
        for r in chr_regions[c]
            if r[4] == 0 || r[5] == 0
                if haskey(ai_losses, "$c")
                    add_seg!(ai_losses["$c"], r[1], r[2], 'a')
                else
                    ai_losses["$c"] = [(r[1], r[2], 1, 0, 0)]
                end
            elseif (r[4] > 1 || r[5] > 1) && r[4] != r[5]
                if haskey(ai_gains, "$c")
                    add_seg!(ai_gains["$c"], r[1], r[2], 'a')
                else
                    ai_gains["$c"] = [(r[1], r[2], 1, 0, 0)]
                end
            end
        end
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

function print_cell(cell::Cell, reporting_data)
    segments = get_segments_for_clean_cell()
    for m in cell.mutation_list
        # Dispatch to chromosomal operation -specific refseq segment generator function, which should alter the provided segment
        # array of ("reference chromosome name", start_bp_of_segment, end_bp_of_segment) tuples according to the mutation operation
        # that was performed. Second parameter contains any data that was returned by the mutation operation.
        m.op(segments, m.data)
    end

    (cumulative_regions, ai_gains, ai_losses, iteration) = reporting_data
    add_segments_to_region_array(segments, cumulative_regions)
    add_segments_to_ai_arrays(segments, ai_gains, ai_losses)
    if haskey(options, "-b")
        write_segments_to_bed_file(segments, iteration)
    end
end


function write_bed_file(chr_regions)
    fname = options["-o"]
    open(fname, "w") do f
        if !haskey(options, "--parent-counts")
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
    println("Wrote Bed file to '$fname'")
end

function get_region_name(chr, position)
    for (pos, loc) in ref_cytobands[chr]
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

function doublecheck_k(cell::Cell)
    k::Int32 = 1
    for gene in driver_genes
        if gene.suppressor == true
            k += 2
        else
            k -= 2
        end
    end
    for chr in cell.chromosomes
        for driver in chr.drivers
            if driver.gene.suppressor == true
                k -= 1
            else
                k += 1
            end
        end
    end
    println("")
    println("--")
    println("Cell with highest fitness (k = $k) evolved through the following chromosomal events:")
end

function print_single_cell(cell::Cell)
    doublecheck_k(cell)
    segments = get_segments_for_clean_cell()
    for m in cell.mutation_list
        # Dispatch to chromosomal operation -specific refseq segment generator function, which should alter the provided segment
        # array of ("reference chromosome name", start_bp_of_segment, end_bp_of_segment) tuples according to the mutation operation
        # that was performed. Second parameter contains any data that was returned by the mutation operation.
        m.op(segments, m.data)
    end
    print_chromosome_regions(segments)
end

function simulate_from_single_cell(pop::Population)
    i = 1
    max_k_idx = 1
    max_k = 1
    cell_count = 1
    mutation_count = 0
    while mutation_count < m
        (pop, max_k_idx) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        if cell_count == 0
            return 0
        end
        max_k = pop.individuals[max_k_idx].k
        mutation_count = size(pop.individuals[max_k_idx].mutation_list)[1]
        print("\rGeneration $i, $cell_count cells, maximum k was $max_k in a cell having $mutation_count mutation events  ")
        i += 1
    end
    print_single_cell(pop.individuals[max_k_idx])
    return cell_count
end

function run_sim_from_single_cell()
    while true
        if simulate_from_single_cell(Population([get_clean_cell(generate_drivers())])) > 0
            println("Done")
            break
        else
            println("")
            println("Crypt died, trying again")
        end
    end
end

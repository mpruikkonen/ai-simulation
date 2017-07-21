#!/home/mipr/julia/julia-903644385b/bin/julia

# you can also run this by include('sim.jl') in the Julia REPL
# and inspect/manipulate the resulting data structures after the
# initial simulation run

# the binary (or source) package for Julia v. 0.6 is required, the
# language syntax and keywords are different in the earlier versions
# that are currently available in e.g. Linux distribution packages

const s = 0.004 # average selective growth advantage of a driver mutation (Bozic/Vogelstein/Nowak 2010, doi:10.1073/pnas.1010978107)
const u = 0.01 # probability of chromosomal event per chromosome per cell division (Lengauer/Kinzler/Vogelstein 1997, doi:10.1038/386623a0)

include("types.jl")
include("whole_chromosome_duplication.jl")
include("chromosome_end_deletion.jl")

# TODO fill in with relevant genes
hg38_driver_genes = [Driver_Gene("TP53", "17", 7661779, 7687550, true),
                     Driver_Gene("SMAD4", "18", 51028394, 51085045, true),
                     Driver_Gene("HNF4A", "20", 44355700, 44434596, false),
                     Driver_Gene("APC", "5", 112707498, 112846239, true),
                     Driver_Gene("KLF5", "13", 73054976, 73077542, false),
                     Driver_Gene("USP12", "13", 27066142, 27171896, false),
                     Driver_Gene("RUNX3", "1", 24899511, 24965157, true),
                     Driver_Gene("TLR3", "4", 186069152, 186088069, true),
                     Driver_Gene("SOX9", "17", 72121020, 72126420, true),
                     Driver_Gene("MYC", "8", 127735434, 127741434, false),
                     Driver_Gene("BMPR1A", "10", 86755786, 86932838, true)]

hg38_chromosome_sizes = Dict("1" => 248956422, "2" => 242193529, "3" => 198295559, "4" => 190214555, "5" => 181538259, "6" => 170805979, "7" => 159345973, "8" => 145138636, "9" => 138394717, "10" => 133797422, "11" => 135086622, "12" => 133275309, "13" => 114364328, "14" => 107043718, "15" => 101991189, "16" => 90338345, "17" => 83257441, "18" => 80373285, "19" => 58617616, "20" => 64444167, "21" => 46709983, "22" => 50818468, "X" => 156040895, "Y" => 57227415)

function chromosome_mutation(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    # Different chromosomal operations go here
    # TODO insertions (keep track which chromosome they are attached to)
    # TODO tandem fusion
    # TODO deletions between two double-strand breaks
    ops = [(0.20, whole_chromosome_duplication),
           (0.60, chromosome_start_deletion),
           (1.00, chromosome_end_deletion)]
    p = rand()
    for (lim, op) in ops
        if(p < lim)
            return op, op(cell, chromosome, push_c1!, push_c2!)
        end
    end
end

function cell_division(cell::Cell)
    new_chromosomes1 = Array{Chromosome, 1}(0)
    new_chromosomes2 = Array{Chromosome, 1}(0)
    new_ml1 = copy(cell.mutation_list)
    new_ml2 = copy(cell.mutation_list)
    mcount1 = mcount2 = 0
    chromosome_idx1 = chromosome_idx2 = 1

    for chromosome in cell.chromosomes
        # Assume we start with CIN, so chromosomal event probability u = 10^-2 per
        # chromosome per cell division according to Lengauer/Kinzler/Vogelstein (1997).
        # The assumption of initial/very early CIN in CRC is reasonable according to Nowak's
        # Evolutionary Dynamics (chapter 12.5), but backed by little experimental evidence
        # and the effective populations in crypts are maybe not as small as Nowak speculated
        # if dedifferentiation plays a role, so we probably want to simulate other initial
        # scenarios as well at some point.
        if(rand() < u)
            cidx1_increase = cidx2_increase = 0
            function push_c1!(chromosome::Chromosome)
                push!(new_chromosomes1, chromosome)
                cidx1_increase += 1
            end
            function push_c2!(chromosome::Chromosome)
                push!(new_chromosomes2, chromosome)
                cidx2_increase += 1
            end

            # Dispatch to operation-specific functions that should call the provided push_c callbacks
            # as many times as the operation produces chromosomes for each child. Return values should
            # contain the number of fitness-increasing driver mutations that occurred for each child
            # (can be negative). Third return value is passed as-is to refseq segment generation.
            # TODO whole genome duplication / multiple-chromosome events?
            (op, (mc1, mc2, data)) = chromosome_mutation(cell, chromosome, push_c1!, push_c2!)

            mcount1 += mc1; mcount2 += mc2
            push!(new_ml1, Mutation(op, chromosome_idx1, true, data))
            push!(new_ml2, Mutation(op, chromosome_idx2, false, data))
            chromosome_idx1 += cidx1_increase
            chromosome_idx2 += cidx2_increase
        else
            push!(new_chromosomes1, chromosome)
            push!(new_chromosomes2, chromosome)
            chromosome_idx1 += 1
            chromosome_idx2 += 1
        end
    end

    return Cell(new_chromosomes1, cell.k+mcount1, new_ml1), Cell(new_chromosomes2, cell.k+mcount2, new_ml2)
end

function get_next_generation(pop::Population)
    a = Array{Cell}(0)
    max_k_idx = 1
    for cell in pop.individuals
        # cell dies with probability 0.5*(1-s)^k (Bozic/Vogelstein/Nowak 2010)
        if(rand() >= (1-s)^cell.k/2)
            (child1, child2) = cell_division(cell)
            push!(a, child1)
            push!(a, child2)
            if(child1.k > a[max_k_idx].k)
                max_k_idx = size(a)[1]-1
            end
            if(child2.k > a[max_k_idx].k)
                max_k_idx = size(a)[1]
            end
        end
    end
    return Population(a), max_k_idx
end

function print_chromosomes(segments)
    println("--")
    println("Final chromosomal layout of the cell with the highest fitness:")
    idx = 1
    for chr in segments
        print("Chromosome $idx:")
        for seg in chr
            print(" [$(seg[1]): $(seg[2])-$(seg[3])]")
        end
        println("")
        idx += 1
    end
end

function print_mutation(idx, old_seg, new_seg)
    print("Chromosome $idx: ")
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
        push!(segments, [("$i", 0, hg38_chromosome_sizes["$i"])])
        push!(segments, [("$i", 0, hg38_chromosome_sizes["$i"])])
    end
    push!(segments, [("X", 0, hg38_chromosome_sizes["X"])])
    push!(segments, [("Y", 0, hg38_chromosome_sizes["Y"])])
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
    
    print_chromosomes(segments)
end

function simulate(pop::Population)
    i = 1
    max_k_idx = 1
    cell_count = 1
    while cell_count > 0 && cell_count < 10000
        println("Generation $i, $cell_count cells, maximum number of driver mutations in a cell: $(pop.individuals[max_k_idx].k)")
        (pop, max_k_idx) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        i += 1
    end
    if cell_count > 0
        println("--")
        println("Cell with highest fitness evolved through the following chromosomal events:")
        print_cell(pop.individuals[max_k_idx])
    end
    return cell_count
end

function get_clean_cell(driver_genes)
    chromosomes = Array{Chromosome, 1}(0)
    for i in 1:22
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["$i"], driver_genes["$i"]))
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["$i"], driver_genes["$i"]))
    end

    push!(chromosomes, Chromosome(hg38_chromosome_sizes["X"], driver_genes["X"]))
    push!(chromosomes, Chromosome(hg38_chromosome_sizes["Y"], driver_genes["Y"]))

    # start with CIN and count it as a driver mutation
    return Cell(chromosomes, 1, [])
end

function generate_drivers()
    d = Dict{String, Array}()
    for i in 1:22
        d["$i"] = []
    end
    d["X"] = []
    d["Y"] = []
    # just take the middle of the gene as the gene position for now (never split genes in chromosome breakages)
    for gene in hg38_driver_genes
        push!(d[gene.native_chromosome], Driver_Locus(floor((gene.native_startpos+gene.native_endpos)/2), gene))
    end
    return d
end

function run_sim_from_single_cell()
    while true
        if simulate(Population([get_clean_cell(generate_drivers())])) > 0
            println("Done")
            break
        else
            println("Crypt died, trying again")
        end
    end
end

function run_sim_for_10m_crypts()
    # TODO? try running with the 10m crypts of the human colon with 1 or 1-4 cell effective populations per crypt
    # - maybe try to model the sequential breakdown of the 1-4 cell limit of the stem cell population as well as the cancer progresses
    # - could simulate just the fraction of the crypts estimated to have CIN or other initial event and introduce new crypts as time/generations progress
    # - could also try to introduce dedifferentiating cells from the differentiated (non-simulated) crypt population back to the simulated cells somehow if data exists on the frequency of such events at different stages of mutation accumulation
    # - set initial conditions for familial adenomatous polyposis and improve simulation by trying to get results to agree with experimental data for that as well?
end

pop = run_sim_from_single_cell()


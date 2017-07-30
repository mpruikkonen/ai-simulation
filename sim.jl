#!/home/mipr/julia/julia-903644385b/bin/julia

# you can also run this by include('sim.jl') in the Julia REPL
# and inspect/manipulate the resulting data structures after the
# initial simulation run

# the binary (or source) package for Julia v. 0.6 is required, the
# language syntax and keywords are different in the earlier versions
# that are currently available in e.g. Linux distribution packages

const s = 0.004 # average selective growth advantage of a driver mutation (Bozic/Vogelstein/Nowak 2010, doi:10.1073/pnas.1010978107)
const u = 0.01 # probability of chromosomal event per chromosome per cell division (Lengauer/Kinzler/Vogelstein 1997, doi:10.1038/386623a0)

options = Dict("-o" => "out.bed")

include("types.jl")
include("genes.jl")
include("chromosome_data.jl")
include("reporting.jl")
for file in filter(x -> endswith(x, ".jl"), readdir("operations"))
    include("operations/$file")
end

function chromosome_mutation(cell::Cell, chromosome::Chromosome, push_c1!::Function, push_c2!::Function)
    # Different chromosomal operations go here
    # TODO whole genome duplication
    # TODO chromosome loss
    # TODO end insert
    ops = [(0.20, whole_chromosome_duplication),
           (0.40, chromosome_tandem_duplication),
           (0.60, chromosome_start_deletion),
           (0.80, chromosome_mid_deletion),
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

function simulate(pop::Population)
    i = 1
    max_k_idx = 1
    cell_count = 1
    while cell_count > 0 && cell_count < 10000
        print("\rGeneration $i, $cell_count cells, maximum number of driver mutations in a cell: $(pop.individuals[max_k_idx].k)")
        (pop, max_k_idx) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        i += 1
    end
    if cell_count > 0
        println("")
        println("--")
        println("Cell with highest fitness evolved through the following chromosomal events:")
        print_cell(pop.individuals[max_k_idx])
    end
    return cell_count
end

function get_centromere_pos(chr)
    ppos = 0
    for (pos, loc) in hg38_cytobands[chr]
        if(loc[1] == 'q')
            return ppos
        end
        ppos = pos
    end
end

function get_clean_cell(driver_genes)
    chromosomes = Array{Chromosome, 1}(0)
    for i in 1:22
        cpos = get_centromere_pos("$i")
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["$i"], cpos, driver_genes["$i"]))
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["$i"], cpos, driver_genes["$i"]))
    end

    push!(chromosomes, Chromosome(hg38_chromosome_sizes["X"], get_centromere_pos("X"), driver_genes["X"]))
    push!(chromosomes, Chromosome(hg38_chromosome_sizes["Y"], get_centromere_pos("Y"), driver_genes["Y"]))

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
            println("")
            println("Crypt died, trying again")
        end
    end
end

function print_help()
    print("""
Usage:
  sim.jl [OPTION]...

Options:
  -h, --help         print this help
  --std-bed          do not output parent-specific counts to BED file
  -o, --output file  name of BED output file (default: '$(options["-o"])')
""")
    quit()
end

function parse_cmdline()
    while !isempty(ARGS)
        a = shift!(ARGS)
        if a == "--std-bed"
            options[a] = ""
        elseif (a == "--output" || a == "-o") && !isempty(ARGS)
            options["-o"] = shift!(ARGS)
        elseif a == "--help" || a == "-h"
            print_help()
        else
            println("Invalid option '$a'")
            print_help()
        end
    end
end

parse_cmdline()
pop = run_sim_from_single_cell()

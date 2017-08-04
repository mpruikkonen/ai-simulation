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
include("chromosome_operations.jl")
include("reporting.jl")

function cell_division(cell::Cell)
    # Assume we start with CIN, so chromosomal event probability u = 10^-2 per
    # chromosome per cell division according to Lengauer/Kinzler/Vogelstein (1997).
    # The assumption of initial/very early CIN in CRC is reasonable according to Nowak's
    # Evolutionary Dynamics (chapter 12.5), but backed by little experimental evidence
    # and the effective populations in crypts are maybe not as small as Nowak speculated
    # if dedifferentiation plays a role, so we probably want to simulate other initial
    # scenarios as well at some point.
    ccount = size(cell.chromosomes)[1]
    if rand() >= (1-u)^ccount
        return chromosome_mutation(cell), cell
    else
        return cell, cell
    end
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
    max_k = 1
    cell_count = 1
    mutation_count = 0
    while cell_count > 0 && mutation_count < 39
        (pop, max_k_idx) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        if cell_count > 0
            max_k = pop.individuals[max_k_idx].k
            mutation_count = size(pop.individuals[max_k_idx].mutation_list)[1]
        end
        print("\rGeneration $i, $cell_count cells, maximum k was $max_k in a cell having $mutation_count mutation events  ")
        i += 1
    end
    if cell_count > 0
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

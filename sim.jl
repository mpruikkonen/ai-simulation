#!/home/mipr/julia/julia-903644385b/bin/julia

# you can also run this by include('sim.jl') in the Julia REPL
# and inspect/manipulate the resulting data structures after the
# initial simulation run

# the binary (or source) package for Julia v. 0.6 is required, the
# language syntax and keywords are different in the earlier versions
# that are currently available in e.g. Linux distribution packages

function print_help()
    print("""
Usage:
  sim.jl [OPTION]...

Options:
  -h, --help            print this help
  -p, --param file      name of parameter file (default: '$(options["-p"])')
  -g, --gains file      name of BedGraph AI gains output file (default: '$(options["-g"])')
  -l, --losses file     name of BedGraph AI losses output file (default: '$(options["-l"])')
  -c, --cumulative file name of BedGraph cumulative output file
  -b, --bed-out dir     in addition to BedGraph files, output individual Bed files for each simulation iteration to dir
  -i, --iterations n    number of simulation iterations (default: '$(options["-i"])')
  -K, --k-iterations n  number of iterations to determine k limit for a simulation iteration (default: '$(options["-K"])')
  -k, --k-limit n       use precalculated k limit n as simulation stopping condition instead of determining it by simulation
  -r, --rng-seed n      use a specific rng seed (to replicate an earlier run)
  -s, --single          run single iteration of simulation from single cell, verbose output (old behavior, useful for debug)
      --parent-counts   output also parent-specific counts to single iteration output file
""")
    quit()
end

options = Dict("-g" => "ai_gains.bedgraph", "-l" => "ai_losses.bedgraph", "-i" => 1000, "-K" => 50, "-r" => Dates.datetime2epochms(now()), "-p" => "ai_curated.jl")

function parse_cmdline()
    while !isempty(ARGS)
        a = shift!(ARGS)
        if a == "--parent-counts"
            options[a] = ""
        elseif (a == "--single" || a == "-s")
            options["-s"] = ""
        elseif (a == "--param" || a == "-p") && !isempty(ARGS)
            options["-p"] = shift!(ARGS)
        elseif (a == "--gains" || a == "-g") && !isempty(ARGS)
            options["-g"] = shift!(ARGS)
        elseif (a == "--losses" || a == "-l") && !isempty(ARGS)
            options["-l"] = shift!(ARGS)
        elseif (a == "--cumulative" || a == "-c") && !isempty(ARGS)
            options["-c"] = shift!(ARGS)
        elseif (a == "--bed-out" || a == "-b") && !isempty(ARGS)
            options["-b"] = shift!(ARGS)
        elseif a == "--iterations" || a == "-i" && !isempty(ARGS)
            options["-i"] = parse(Int64, shift!(ARGS))
        elseif a == "--k-iterations" || a == "-K" && !isempty(ARGS)
            options["-K"] = parse(Int64, shift!(ARGS))
        elseif a == "--k-limit" || a == "-k" && !isempty(ARGS)
            options["-k"] = parse(Float64, shift!(ARGS))
        elseif a == "--rng-seed" || a == "-r" && !isempty(ARGS)
            options["-r"] = parse(Int64, shift!(ARGS))
        elseif a == "--help" || a == "-h"
            print_help()
        else
            println("Invalid option '$a'")
            print_help()
        end
    end
    rngseed = options["-r"]
    srand(rngseed)
    println("rng seed is $rngseed")
end

parse_cmdline()

include("types.jl")
include("util.jl")
for file in filter(x -> endswith(x, ".jl"), readdir("operations"))
    include("operations/$file")
end
include("""parameters/$(options["-p"])""")
include("reporting.jl")

const op_psum = check_op_psum(ops)

function chromosome_mutation(cell::Cell)
    r = rand() * op_psum
    for (p, op) in ops
        if(r < p)
            (c, k, data) = op(cell)
            ml = copy(cell.mutation_list)
            push!(ml, Mutation(op, data))
            return Cell(c, k, ml)
        end
        r -= p
    end
end

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

function simulate_klimit(pop::Population)
    i = 1
    max_k_idx = 1
    max_k = 1
    cell_count = 1
    mutation_count = 0
    while mutation_count < m
        (pop, max_k_idx) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        if cell_count == 0
            return 0, 0
        end
        max_k = pop.individuals[max_k_idx].k
        mutation_count = size(pop.individuals[max_k_idx].mutation_list)[1]
        i += 1
    end
    return cell_count, max_k
end

function simulate(pop::Population, k_limit, reporting_data)
    i = 1
    max_k_idx = 1
    max_k = 1
    cell_count = 1
    while max_k < k_limit
        (pop, max_k_idx) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        if cell_count == 0
            return 0
        end
        max_k = pop.individuals[max_k_idx].k
        i += 1
    end
    print_cell(pop.individuals[max_k_idx], reporting_data)
    return cell_count
end

function get_centromere_pos(chr)
    ppos = 0
    for (pos, loc) in ref_cytobands[chr]
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
        push!(chromosomes, Chromosome(ref_chromosome_sizes["$i"], cpos, driver_genes["$i"]))
        push!(chromosomes, Chromosome(ref_chromosome_sizes["$i"], cpos, driver_genes["$i"]))
    end

    push!(chromosomes, Chromosome(ref_chromosome_sizes["X"], get_centromere_pos("X"), driver_genes["X"]))
    push!(chromosomes, Chromosome(ref_chromosome_sizes["Y"], get_centromere_pos("Y"), driver_genes["Y"]))

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
    for gene in driver_genes
        push!(d[gene.native_chromosome], Driver_Locus(floor((gene.native_startpos+gene.native_endpos)/2), gene))
    end
    return d
end

function get_avg_k(drivers)
    i = 1
    avg_k::Float64 = 0
    iterations = options["-K"]
    println("Simulating for $iterations iterations to get average k for maximally fit cell with at least $m mutation operations")
    while i < iterations
        (cell_count, k) = simulate_klimit(Population([get_clean_cell(drivers)]))
        if cell_count > 0
            print("$k ")
            avg_k += k
            i += 1
        end
    end
    avg_k /= iterations
    println("/ $iterations = $avg_k")
    return avg_k
end

function run_sim()
    cumulative_regions = Dict{String, Array{Tuple{Int, Int, Int, Int, Int}}}()
    ai_gains = Dict{String, Array{Tuple{Int, Int, Int, Int, Int}}}()
    ai_losses = Dict{String, Array{Tuple{Int, Int, Int, Int, Int}}}()
    drivers = generate_drivers()
    avg_k = haskey(options, "-k") ? options["-k"] : get_avg_k(drivers)
    iterations = options["-i"]
    i = 1
    while i < iterations
        if simulate(Population([get_clean_cell(drivers)]), avg_k, (cumulative_regions, ai_gains, ai_losses, i)) > 0
            print("\rIteration $i/$iterations")
            i += 1
        end
    end
    println("\rSimulated $iterations iterations up to a cell with fitness at least k = $avg_k")
    write_bedgraph_file(options["-g"], ai_gains)
    write_bedgraph_file(options["-l"], ai_losses)
    println("""Wrote allelic imbalance gains and losses to '$(options["-g"])' and '$(options["-l"])'""")
    if haskey(options, "-c")
        write_bedgraph_file(options["-c"], cumulative_regions)
        println("""Wrote cumulative BedGraph file to '$(options["-c"])'""")
    end
    if haskey(options, "-b")
        println("""Bed files for each simulation iteration written to '$(options["-b"])'""")
    end
end

if haskey(options, "-s")
    if endswith(options["-o"], ".bedgraph")
        options["-o"] = options["-o"][1:end-5]
    end
    include("debug.jl")
    run_sim_from_single_cell()
else
    if haskey(options, "-b")
        mkpath(options["-b"])
    end
    run_sim()
end

#!/home/mipr/julia/julia-903644385b/bin/julia

# you can also run this by include('sim.jl') in the Julia REPL
# and inspect/manipulate the resulting data structures after the
# initial simulation run

# the binary (or source) package for Julia v. 0.6 is required, the
# language syntax and keywords are different in the earlier versions
# that are currently available in e.g. Linux distribution packages

struct Gene
    name::String
    native_chromosome::String
    native_startpos::UInt32
    native_endpos::UInt32
end

struct Driver_locus
    pos::UInt32
    s::Real
    gene::Gene
end

struct Chromosome
    size::UInt32
    drivers::Array{Driver_locus, 1}
end

struct Cell
    chromosomes::Array{Chromosome, 1}
    mutation_count::UInt8
    # TODO mutation_list/parent (if we want to actually trace the operations that occurred in particular cells)
end

struct Population
    individuals::Array{Cell, 1}
end

# TODO fill in with relevant genes
# TODO actually use the driver gene data in the simulation
hg38_driver_genes = [Gene("MYC", "8", 127735434, 127741434),
                     Gene("APC", "5", 112707498, 112846239)]

hg38_chromosome_sizes = Dict("1" => 248956422, "2" => 242193529, "3" => 198295559, "4" => 190214555, "5" => 181538259, "6" => 170805979, "7" => 159345973, "8" => 145138636, "9" => 138394717, "10" => 133797422, "11" => 135086622, "12" => 133275309, "13" => 114364328, "14" => 107043718, "15" => 101991189, "16" => 90338345, "17" => 83257441, "18" => 80373285, "19" => 58617616, "20" => 64444167, "21" => 46709983, "22" => 50818468, "X" => 156040895, "Y" => 57227415)

function get_next_generation(pop::Population)
    a = Array{Cell}(0)
    max_mut = 0
    for cell in pop.individuals
        new_chromosomes1 = Array{Chromosome, 1}(0)
        new_chromosomes2 = Array{Chromosome, 1}(0)
        mcount = 0
        for chromosome in cell.chromosomes
            # Assume we start with CIN, so chromosomal event probability = 10^-2 per
            # chromosome per cell division according to Lengauer/Kinzler/Vogelstein (1997).
            # The assumption of initial/very early CIN in CRC is reasonable according to Nowak's
            # Evolutionary Dynamics (chapter 12.5), but backed by little experimental evidence
            # and the effective populations in crypts are maybe not as small as Nowark speculated
            # if dedifferentiation plays a role, so we probably want to simulate other initial
            # scenarios as well at some point.
            if(rand() < 0.01)
                mcount += 1
                # different chromosomal operations go here, also possible to
                # include missegregations etc. that affect both children
                # TODO insertions (keep track which chromosome they are attached to)
                # TODO tandem fusion
                if(rand() < 0.25)
                    #println("Whole-chromosome duplication")
                    push!(new_chromosomes1, chromosome)
                    push!(new_chromosomes1, chromosome)
                else
                    # TODO keep track of centromere positions and check which half survives?
                    # TODO deletions between two double-strand breaks
                    delete_bases::UInt32 = floor(rand()*chromosome.size)
                    new_drivers = Array{Driver_locus, 1}(0)
                    if(rand() < 0.5)
                        #println("Deletion of $delete_bases/$(chromosome.size) bases at beginning of chromosome")
                        for driver in chromosome.drivers
                            if(driver.pos > delete_bases)
                                push!(new_drivers, Driver_locus(driver.pos-delete_bases, driver.s, driver.gene))
                            end
                        end
                        push!(new_chromosomes1, Chromosome(chromosome.size-delete_bases, new_drivers))
                    else
                        #println("Deletion of $delete_bases/$(chromosome.size) bases at end of chromosome")
                        for driver in chromosome.drivers
                            if(driver.pos < delete_bases)
                                push!(new_drivers, Driver_locus(driver.pos, driver.s, driver.gene))
                            end
                        end
                        push!(new_chromosomes1, Chromosome(chromosome.size-delete_bases, new_drivers))
                    end
                end
            else
                push!(new_chromosomes1, chromosome)
            end
            push!(new_chromosomes2, chromosome)
        end
        push!(a, Cell(new_chromosomes1, cell.mutation_count+mcount))
        push!(a, Cell(new_chromosomes2, cell.mutation_count))
        if(cell.mutation_count+mcount > max_mut)
            max_mut = cell.mutation_count+mcount
        end
    end
    return Population(a), max_mut
end

function simulate(pop)
    i = 1
    max_mut = 0
    while true
        println("Generation $i, $(size(pop.individuals)[1]) cells, maximum number of chromosomal mutations in a cell: $max_mut")
        (pop, max_mut) = get_next_generation(pop)
        if(max_mut >= 15)
            return pop
        end
        i += 1
    end
end

function get_clean_cell(driver_genes)
    chromosomes = Array{Chromosome, 1}(0)
    for i in 1:22
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["$i"], driver_genes["$i"]))
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["$i"], driver_genes["$i"]))
    end

    push!(chromosomes, Chromosome(hg38_chromosome_sizes["X"], driver_genes["X"]))
    if(rand() < 0.5)
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["X"], driver_genes["X"]))
    else
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["Y"], driver_genes["Y"]))
    end
    
    return Cell(chromosomes, 0)
end

function generate_drivers()
    d = Dict{String, Array}()
    for i in 1:22
        d["$i"] = []
    end
    d["X"] = []
    d["Y"] = []
    # just take the middle of the gene as the gene position for now (never split genes in chromosome breakages)
    # and use a constant fitness advantage of 0.4% for all driver mutations (Bozic/Vogelstein/Nowak 2010)
    for gene in hg38_driver_genes
        push!(d[gene.native_chromosome], Driver_locus(floor((gene.native_startpos+gene.native_endpos)/2), 1.004, gene))
    end
    return d
end

function run_exponential_sim_from_single_cell()
    simulate(Population([get_clean_cell(generate_drivers())]))
end

function run_sim_for_10m_crypts()
    # TODO? try running with the 10m crypts of the human colon with 1 or 1-4 cell effective populations per crypt
    # - maybe try to model the sequential breakdown of the 1-4 cell limit of the stem cell population as well as the cancer progresses
    # - could simulate just the fraction of the crypts estimated to have CIN or other initial event and introduce new crypts as time/generations progress
    # - could also try to introduce dedifferentiating cells from the differentiated (non-simulated) crypt population back to the simulated cells somehow if data exists on the frequency of such events at different stages of mutation accumulation
    # - set initial conditions for familial adenomatous polyposis and improve simulation by trying to get results to agree with experimental data for that as well?
end

pop = run_exponential_sim_from_single_cell()


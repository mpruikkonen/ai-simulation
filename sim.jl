#!/home/mipr/julia/julia-903644385b/bin/julia

# you can also run this by include('sim.jl') in the Julia REPL
# and inspect/manipulate the resulting data structures after the
# initial simulation run

# the binary (or source) package for Julia v. 0.6 is required, the
# language syntax and keywords are different in the earlier versions
# that are currently available in e.g. Linux distribution packages

const s = 0.004 # average selective growth advantage of a driver mutation (Bozic/Vogelstein/Nowak 2010, doi:10.1073/pnas.1010978107)
const u = 0.01 # probability of chromosomal event per chromosome per cell division (Lengauer/Kinzler/Vogelstein 1997, doi:10.1038/386623a0)

struct Driver_Gene
    name::String
    native_chromosome::String
    native_startpos::UInt32
    native_endpos::UInt32
    suppressor::Bool
end

struct Driver_Locus
    pos::UInt32
    gene::Driver_Gene
end

struct Chromosome
    size::UInt32
    drivers::Array{Driver_Locus, 1}
end

struct Cell
    chromosomes::Array{Chromosome, 1}
    k::UInt8 # number of driver mutations, fitness = (1+s)^k
    # TODO mutation_list/parent (if we want to actually trace the operations that occurred in particular cells)
end

struct Population
    individuals::Array{Cell, 1}
end

# TODO fill in with relevant genes
hg38_driver_genes = [Driver_Gene("MYC", "8", 127735434, 127741434, false),
                     Driver_Gene("APC", "5", 112707498, 112846239, true)]

hg38_chromosome_sizes = Dict("1" => 248956422, "2" => 242193529, "3" => 198295559, "4" => 190214555, "5" => 181538259, "6" => 170805979, "7" => 159345973, "8" => 145138636, "9" => 138394717, "10" => 133797422, "11" => 135086622, "12" => 133275309, "13" => 114364328, "14" => 107043718, "15" => 101991189, "16" => 90338345, "17" => 83257441, "18" => 80373285, "19" => 58617616, "20" => 64444167, "21" => 46709983, "22" => 50818468, "X" => 156040895, "Y" => 57227415)

function whole_chromosome_duplication(cell::Cell, chromosome::Chromosome, new_chromosomes1::Array{Chromosome, 1}, new_chromosomes2::Array{Chromosome, 1})
    #println("Whole-chromosome duplication")
    push!(new_chromosomes1, chromosome)
    push!(new_chromosomes1, chromosome)
    push!(new_chromosomes2, chromosome)
    mcount = 0
    for driver in chromosome.drivers
        if driver.gene.suppressor == false
            mcount += 1 # TODO do we want to count just the first duplication of an oncogene instead of each one?
        end
    end
    return mcount, 0
end

function get_deletion_mcount(cell::Cell, gene::Driver_Gene)
    gcount = 0
    for chromosome in cell.chromosomes
        for driver in chromosome.drivers
            if driver.gene == gene
                gcount += 1
            end
        end
    end

    if gene.suppressor == true
        return gcount == 1 ? 1 : 0 # fitness increases when losing last copy of suppressor
    else
        return gcount > 2 ? -1 : 0 # fitness decreases when losing extra copies of oncogene
    end
end

function chromosome_end_deletion(cell::Cell, chromosome::Chromosome, new_chromosomes1::Array{Chromosome, 1}, new_chromosomes2::Array{Chromosome, 1})
    # TODO keep track of centromere positions and check which half survives?
    # TODO deletions between two double-strand breaks
    delete_bases::UInt32 = floor(rand()*chromosome.size)
    new_drivers = Array{Driver_Locus, 1}(0)
    mcount = 0
    if(rand() < 0.5)
        #println("Deletion of $delete_bases/$(chromosome.size) bases at beginning of chromosome")
        for driver in chromosome.drivers
            if(driver.pos > delete_bases)
                push!(new_drivers, Driver_Locus(driver.pos-delete_bases, driver.gene))
            else
                mcount += get_deletion_mcount(cell, driver.gene)
            end
        end
        push!(new_chromosomes1, Chromosome(chromosome.size-delete_bases, new_drivers))
    else
        #println("Deletion of $delete_bases/$(chromosome.size) bases at end of chromosome")
        for driver in chromosome.drivers
            if(driver.pos < delete_bases)
                push!(new_drivers, Driver_Locus(driver.pos, driver.gene))
            else
                mcount += get_deletion_mcount(cell, driver.gene)
            end
        end
        push!(new_chromosomes1, Chromosome(chromosome.size-delete_bases, new_drivers))
    end
    push!(new_chromosomes2, chromosome)
    return mcount, 0
end

function chromosome_mutation(cell::Cell, chromosome::Chromosome, new_chromosomes1::Array{Chromosome, 1}, new_chromosomes2::Array{Chromosome, 1})
    # different chromosomal operations go here, also possible to include missegregations etc. that affect both children
    # TODO insertions (keep track which chromosome they are attached to)
    # TODO tandem fusion
    ops = [(0.25, whole_chromosome_duplication),
           (1.00, chromosome_end_deletion)]
    p = rand()
    for (lim, op) in ops
        if(p < lim)
            return op(cell, chromosome, new_chromosomes1, new_chromosomes2)
        end
    end
end

function cell_division(cell::Cell)
    # TODO whole genome duplication / multiple-chromosome events?
    new_chromosomes1 = Array{Chromosome, 1}(0)
    new_chromosomes2 = Array{Chromosome, 1}(0)
    mcount1 = 0
    mcount2 = 0
    for chromosome in cell.chromosomes
        # Assume we start with CIN, so chromosomal event probability = 10^-2 per
        # chromosome per cell division according to Lengauer/Kinzler/Vogelstein (1997).
        # The assumption of initial/very early CIN in CRC is reasonable according to Nowak's
        # Evolutionary Dynamics (chapter 12.5), but backed by little experimental evidence
        # and the effective populations in crypts are maybe not as small as Nowak speculated
        # if dedifferentiation plays a role, so we probably want to simulate other initial
        # scenarios as well at some point.
        if(rand() < u)
            (mc1, mc2) = chromosome_mutation(cell, chromosome, new_chromosomes1, new_chromosomes2)
            mcount1 += mc1
            mcount2 += mc2
        else
            push!(new_chromosomes1, chromosome)
            push!(new_chromosomes2, chromosome)
        end
    end
    return Cell(new_chromosomes1, cell.k+mcount1), Cell(new_chromosomes2, cell.k+mcount2), mcount1 > mcount2 ? mcount1 : mcount2
end

function get_next_generation(pop::Population)
    a = Array{Cell}(0)
    max_k = 0
    for cell in pop.individuals
        # cell dies with probability 0.5*(1-s)^k (Bozic/Vogelstein/Nowak 2010)
        if(rand() >= (1-s)^cell.k/2)
            (child1, child2, mcount) = cell_division(cell)
            push!(a, child1)
            push!(a, child2)
            if(cell.k+mcount > max_k)
                max_k = cell.k+mcount
            end
        end
    end
    return Population(a), max_k
end

function simulate(pop)
    i = 1
    max_k = 0
    cell_count = 1
    while cell_count > 0 && cell_count < 10000
        println("Generation $i, $cell_count cells, maximum number of driver mutations in a cell: $max_k")
        (pop, max_k) = get_next_generation(pop)
        cell_count = size(pop.individuals)[1]
        i += 1
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
    if(rand() < 0.5)
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["X"], driver_genes["X"]))
    else
        push!(chromosomes, Chromosome(hg38_chromosome_sizes["Y"], driver_genes["Y"]))
    end

    # start with CIN and count it as a driver mutation
    return Cell(chromosomes, 1)
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


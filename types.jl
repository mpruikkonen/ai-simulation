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
    centromere_pos::UInt32
    drivers::Array{Driver_Locus, 1}
end

struct Mutation
    op::Function
    data
end

struct Cell
    chromosomes::Array{Chromosome, 1}
    k::Int32 # number of driver mutations, fitness = (1+s)^k
    mutation_list::Array{Mutation, 1}
end

struct Population
    individuals::Array{Cell, 1}
end

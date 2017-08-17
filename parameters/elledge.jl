include("ref/hg19.jl")
include("ops/noarmlength.jl")
include("genes/elledge.jl")

const s = 0.002 # average selective growth advantage of a driver mutation (Bozic/Vogelstein/Nowak 2010, doi:10.1073/pnas.1010978107)
const u = 0.01 # probability of chromosomal event per chromosome per cell division (Lengauer/Kinzler/Vogelstein 1997, doi:10.1038/386623a0)
const m = 39 # average number of mutation operations per simulation endpoint (Zack/Beroukhim 2013, doi:10.1038/ng.2760)

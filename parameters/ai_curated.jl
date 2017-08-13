const s = 0.004 # average selective growth advantage of a driver mutation (Bozic/Vogelstein/Nowak 2010, doi:10.1073/pnas.1010978107)
const u = 0.01 # probability of chromosomal event per chromosome per cell division (Lengauer/Kinzler/Vogelstein 1997, doi:10.1038/386623a0)
const m = 39 # average number of mutation operations per simulation endpoint (Zack/Beroukhim 2013, doi:10.1038/ng.2760)

include("ref/hg19.jl")

# 12% of focal amplifications and 26% of focal deletions are telomere-bound according to Zack/Beroukhim 2013 (doi:10.1038/ng.2760),
# whereas focal amplifications and deletions are relatively balanced (median of 11 amplifications and 12 deletions per sample).
# Frequency for whole genome duplications in CRC is stated as .43/sample. Following approximates these relative frequencies.

# Zack/Beroukhim also mentions arm-length or longer SCNAs (median 3 amplifications and 5 deletions per sample) and copy-neutral
# LOH events (median 1 per sample) (TODO: should these be implemented separately as well?)
# For now just estimate frequency of whole-chromosome events to amount to .57 per sample with 3 duplications for 5 deletions.

const ops = [(.43 * 1/24, whole_genome_duplication),            # 0.018
             (3/8 * .57 * 1/24, whole_chromosome_duplication),  # 0.009
             (5/8 * .57 * 1/24, whole_chromosome_loss),         # 0.015
             (.88 * 11/24, chromosome_tandem_duplication),      # 0.403
             (.06 * 11/24, chromosome_start_insertion),         # 0.028
             (.06 * 11/24, chromosome_end_insertion),           # 0.028
             (.74 * 12/24, chromosome_mid_deletion),            # 0.370
             (.13 * 12/24, chromosome_start_deletion),          # 0.065
             (.13 * 12/24, chromosome_end_deletion)]            # 0.065

include("ai_curated_genes.jl") # array of curated driver genes from ai paper

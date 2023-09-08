#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
sampleNo <- as.numeric(args[1])

if (sampleNo == '' | is.na(sampleNo)) {
    cat("ERROR: No sample number provided as an argument!\n")
    quit(status=1)
}

cat(paste0('INFO: SampleNo is ', sampleNo, '\n'))

all_peaks_bed_file_list <- 'all_peak_bed_files.txt'
all_files <- readLines(all_peaks_bed_file_list)
samplefile <- all_files[sampleNo]
fragmentfile <- paste0(strsplit(samplefile, split='atac_peaks.bed')[[1]], 'atac_fragments.tsv.gz')





all_chrs <- readLines('all_contigs.txt')



sample <- strsplit(samplefile, split='Multiome/')
sample <- strsplit(sample[[1]], split='/outs')[[2]][1]

print(sampleNo)

print(samplefile)

print(sample)

sample_peaks <- fread(fragmentfile, header=F)


setnames(sample_peaks, c('chr','start','end','id','N'))
setkey(sample_peaks, chr, start, end)

merged_peaks <- fread('merged_peaks.tsv')
setnames(merged_peaks, c('chr','start','end'))
setkey(merged_peaks, chr, start, end)

ov <- foverlaps(sample_peaks, merged_peaks, type='within')
# exclude non-overlapping rows
ov <- ov[!is.na(start)]
ov <- ov[!is.na(end)]

ov <- ov[, list('N'=sum(N)), by=list(chr, start, end, id)]
ov[, chr := factor(chr, levels=all_chrs)]
setkey(ov, chr, start, end)



outfile <- paste0('reassigned_fragments/', sample, '.reassigned_fragments.tsv')
bgzipfile <- paste0('reassigned_fragments/', sample, '.reassigned_fragments.tsv.gz')

# Write header

write(paste0('# id=', sample) ,file=outfile)
write(paste0('# description=', sample) ,file=outfile, append=T)
write('#\n# reference_path=/vf/db/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A', file=outfile, append=T)
write('# reference_fasta_hash=b6f131840f9f337e7b858c3d1e89d7ce0321b243\n# reference_gtf_hash=3b4c36ca3bade222a5b53394e8c07a18db7ebb11', file=outfile, append=T)
write('# reference_version=2020-A\n# mkref_version=cellranger-arc-1.0.0\n#', file=outfile, append=T)
write(paste0('# primary_contig=', all_chrs, collapse='\n'), file=outfile, append=T)

# Write sample
fwrite(ov, file=outfile, quote=F, row.names=F, col.names=F, sep='\t', append=T)

# bgzip
system(paste0('bgzip ', outfile))

# tabix index
system(paste0('tabix -p bed ', bgzipfile))
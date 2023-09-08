#!/usr/bin/env Rscript

if(file.exists('merged_peaks.tsv')) {
    cat("WARNING: merged_peaks.tsv already exists. This script would overwrite it.\n")
    cat("         Remove merged_peaks.tsv and try again if you really want this to happen.\n")
    quit(status=1)
}

library(data.table)
library(foreach)
library(GenomicRanges)

all_files <- readLines('all_peak_bed_files.txt')
all_chrs <- readLines('all_contigs.txt')




for(chr in all_chrs) {
    o <- foreach(file=all_files, .combine='rbind') %do% {
        sample <- strsplit(file, split='Multiome/')
        sample <- strsplit(sample[[1]], split='/outs')[[2]][1]
        dat <- fread(file, header=F)
        setnames(dat, c('CHR','start','end'))
        dat <- dat[CHR==chr]
        dat[, 'sample' := sample]
        return(dat)
    }
    fwrite(o, file=paste0(chr, '_peaks.tsv'), sep='\t', quote=F, row.names=F, col.names=T)
}

# dat.unmerged <- foreach(file=list.files(pattern='*peaks.tsv'), .combine='rbind') %do% {
#     fread(file)
# }
# dat.unmerged[, width := end-start]

merged_peaks_files <- list.files(pattern='*peaks.tsv')

dat.merged <- foreach(file=merged_peaks_files, .combine='rbind') %do% {
    dat <- fread(file)
    dat <- dat[, as.data.table(reduce(IRanges(start, end))), by=list(CHR)]
    return(dat)
}

# dat.unmerged[, sample := NULL]
# dat.unmerged[, type := 'unmerged']
# dat.merged[, type := 'merged']
# dat.all <- rbindlist(list(dat.merged, dat.unmerged))

# Write all-chr peaks file
fwrite(dat.merged[, c('CHR','start','end')], file='merged_peaks.tsv', quote=F, col.names=F, row.names=F, sep='\t')

# remove individual peaks files
file.remove(merged_peaks_files)
quit(status=0)












ggplot(dat.all, aes(x=width, color=type)) +geom_density() + facet_wrap(~CHR)

ggplot(dat.all, aes(x=width, color=type)) +geom_density() + facet_wrap(~type, scales='free')






# Generate keyed data.table with chromosome-specific ranges
bed <- fread('filename.bed')
setnames(bed, c("CHROM","start","stop"))
setkey(bed, CHROM, start, stop)

# To use vcf sites as a "range" create two columns for
# "start" and "stop" that are just the POS value, in reality a 1-length range
vcf <- fread('filename.vcf', header=T, skip="#CHROM")
setnames(vcf, "#CHROM", "CHROM")
vcf[, "stop" := POS]
setnames(vcf, "POS", "start")

# Give each row a unique ID for match (or inverse match)
vcf[, id := 1:.N]

# Perform overlap search, returning vector of ids
# corresponding to sites that are contained within at least one row of bed file,
# e.g. bed_start <= vcf_POS <= bed_stop
ids_within_bed_ranges <- unique(foverlaps(vcf, bed, type="within")[!is.na(start), id])

# Perform subset of interest
vcf.outside <- vcf[! id %in% ids_within_bed_ranges]
vcf.inside <- vcf[id %in% ids_within_bed_ranges]

# Fix vcf header names and remove extra columns
setnames(vcf.outside, "start", "POS")
vcf.outside[, c("stop","id") := NULL]

fn <- 'within.tsv'

sample <- fread('within.tsv')
setnames(sample, c('chr','start','end','id','N'))
setkey(sample, chr, start, end)

consensus <- fread('consensus_merged_peaks.tsv')
setnames(consensus, c('chr','start','end'))
setkey(consensus, chr, start, end)

ov <- foverlaps(sample, consensus, type='within')
# exclude non-overlapping rows
ov <- ov[!is.na(start)]
ov <- ov[!is.na(end)]

ov <- ov[, list('N'=sum(N)), by=list(chr, start, end, id)]
ov[, chr := factor(chr, levels=all_chrs)]
setkey(ov, chr, start, end)


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sequenza))
suppressPackageStartupMessages(library(falcon))

option_list <- list(
    # Required
    make_option(c("--sample_name", "-s"), dest="sample_name", action="store", default = NA, type = 'character',
                help = "[Required] Name of the sample"),
    make_option(c("--sample_input", "-i"), dest="sample_input", action="store", default = NA, type = 'character',
                help = "[Required] Path to the sample.seqz.gz file"),
    make_option(c("--centromeres", "-c"), dest="centromere", action="store", default = NA, type = 'character',
                help = "[Required] Path to file defining centromere locations"),
    make_option(c("--out_dir", "-o"), dest="out_dir", action="store", default = NA, type = 'character',
                help = "[Required] Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

#reading in data
sample_name <- opt[["sample_name"]]
sample_input <- opt[["sample_input"]]
centromere_input <- opt[["centromere"]]
out_dir <- opt[["out_dir"]]

centromeres<-read.csv(centromere_input, sep = '\t', header = FALSE, col.names = c("chromosome","start","stop"))
chromosome <- centromeres$chromosome

#function created using falcon to use the centromere position given
falcon.seg.seqz <- function(data.file, chromosome, centromeres){
    require(sequenza)
    seqz <- read.seqz(data.file, chr_name = chromosome)
    chrom <- gsub(chromosome, pattern = "chr", replacement = "")
    centromeres$chromosome <- gsub(centromeres$chromosome, pattern = "chr", replacement = "")
    seqz   <- seqz[seqz$zygosity.normal == "het", ]
    get.tauhat <- function(seqz, ...) {
        require(falcon)
        at     <- round(seqz$depth.tumor * seqz$Af, 0)
        bt     <- round(seqz$depth.tumor * seqz$Bf, 0)
        an     <- round(seqz$depth.normal * 0.55, 0)
        bn     <- round(seqz$depth.normal * 0.45, 0)
        getChangepoints(data.frame(AT = at, BT = bt, AN = an, BN = bn) , ...)
    }
    p <- seqz$position < centromeres$start.pos[centromeres$chromosome == chrom]
    q <- seqz$position > centromeres$end.pos[centromeres$chromosome == chrom]
    pos.p    <- seqz$position[p]
    pos.q    <- seqz$position[q]
    l.p <- length(pos.p)
    l.q <- length(pos.q)
    do.breaks <- function(chrom, tauhat) {
        start.pos <- tauhat
        start.pos[-1] <- tauhat[-1]+1
        data.frame(chrom = chrom,
                   start.pos = start.pos[-(length(start.pos))],
                   end.pos = tauhat[-1])
    }
    chrom  <- unique(seqz$chromosome)
    if (l.p > 1 & l.q > 1) {
        tauhat.p <- get.tauhat(seqz[p, ], verbose = FALSE)
        tauhat.p <- c(min(pos.p), pos.p[tauhat.p], max(pos.p))
        tauhat.q <- get.tauhat(seqz[q, ], verbose = FALSE)
        tauhat.q <- c(min(pos.q), pos.q[tauhat.q], max(pos.q))
        breaks.p <- do.breaks(chrom, tauhat.p)
        breaks.q <- do.breaks(chrom, tauhat.q)
        rbind(breaks.p, breaks.q)
    } else if (l.p < 2) {
        tauhat.q <- get.tauhat(seqz[q, ], verbose = FALSE)
        tauhat.q <- c(min(pos.q), pos.q[tauhat.q], max(pos.q))
        do.breaks(chrom, tauhat.q)
    } else if (l.q < 2 ) {
        tauhat.p <- get.tauhat(seqz[p, ], verbose = FALSE)
        tauhat.p <- c(min(pos.p), pos.p[tauhat.p], max(pos.p))
        do.breaks(chrom, tauhat.p)
    } else {
        stop("Segmentation went wrong...")
    }
}

seqz.data <- read.seqz(sample_input)
str(seqz.data, vec.len=2)

brk<-list()
for (i in chromosome){
    brk[[i]] <- falcon.seg.seqz(sample_input, chromosome = i, centromeres = centromeres)
}
breaks <- do.call(rbind, brk)

options("scipen"=100, "digits"=4)
#this will force the program to read the regions as only straight values and not as exponential

sample <- sequenza.extract(sample_input, breaks=breaks, chromosome.list=chromosome, parallel = 4)

#estimates cellularity and ploidy
CP <- sequenza.fit(sample, mc.core = 4, segment.filter= 5e6)

sequenza.results(sample, CP, out.dir=out_dir, sample.id = sample_name, chromosome.list = chromosome)

png(file=paste(paste(out_dir,sample_name,sep=''),'Maximum_Likelihood_plot.png',sep=''),type='cairo')
cp.plot(CP)
cp.plot.contours(CP, add = TRUE,
   likThresh = c(0.999, 0.95),
   col = c("lightsalmon", "red"), pch = 20)

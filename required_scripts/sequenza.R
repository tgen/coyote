
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sequenza))
suppressPackageStartupMessages(library(falcon))

option_list <- list(
    # Required
    make_option(c("--sample_name", "-s"), dest="sample_name", action="store", default = NA, type = 'character',
                help = "[Required] Name of the sample"),
    make_option(c("--sample_input", "-i"), dest="sample_input", action="store", default = NA, type = 'character',
                help = "[Required] Path to the sample.seqz.gz file"),
    make_option(c("--out_dir", "-o"), dest="out_dir", action="store", default = NA, type = 'character',
                help = "[Required] Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

#reading in data
sample_name <- opt[["sample_name"]]
sample_input <- opt[["sample_input"]]
out_dir <- opt[["out_dir"]]

#define a dataframe with centromere start and end
chromosome <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chr32", "chr33", "chr34", "chr35", "chr36", "chr37", "chr38", "chrX")

start.pos <- c(69540000,57400000,31240000,35185000,30100000,40470000,78545000,750000,61045000,34025000,26795000,20000,63130000,12025000,42685000,21210000,51325000,25210000,0,45385000,50000,30320000,50740000,5000,19050000,26435000,45825000,41175000,24320000,37875000,35000,35000,31365000,12100000,26455000,0,140000,21640000,49595000)

end.pos <- c(69655000,57940000,31640000,35235000,30115000,40605000,78570000,775000,61055000,34075000,26855000,25000,63235000,12035000,42720000,21250000,51365000,25635000,1000,45390000,60000,30325000,50745000,60000,19085000,26440000,45840000,41180000,24325000,37880000,115000,50000,31370000,12505000,26465000,1000,145000,21645000,49735000)

centromeres<-data.frame(chromosome=chromosome,start.pos=start.pos,end.pos=end.pos)

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
chromosome.list <- c(1:38,"X")
for (i in chromosome.list){
    brk[[i]] <- falcon.seg.seqz(sample_input, chromosome = i, centromeres = centromeres)
}
breaks <- do.call(rbind, brk)

options("scipen"=100, "digits"=4)
#this will force the program to read the regions as only straight values and not as exponential

sample <- sequenza.extract(sample_input, breaks=breaks, chromosome.list=chromosome.list, parallel = 4)

#estimates cellularity and ploidy
CP <- sequenza.fit(sample, mc.core = 4, segment.filter= 5e6)

sequenza.results(sample, CP, out.dir=out_dir, sample.id = sample_name, chromosome.list = chromosome.list)

png(file=paste(paste(out_dir,sample_name,sep=''),'Maximum_Likelihood_plot.png',sep=''),type='cairo')
cp.plot(CP)
cp.plot.contours(CP, add = TRUE,
   likThresh = c(0.999, 0.95),
   col = c("lightsalmon", "red"), pch = 20)


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sequenza))
suppressPackageStartupMessages(library(falcon))

option_list <- list(
    # Required
    make_option(c("--sample_name", "-s"), dest="sample_name", action="store", default = NA, type = 'character',
                help = "[Required] The title used for each plot"),
    make_option(c("--sample_input", "-i"), dest="sample_input", action="store", default = NA, type = 'character',
                help = "[Required] The title used for each plot")
)

opt <- parse_args(OptionParser(option_list=option_list))

#reading in data
sample_name <- opt[["sample_name"]]
sample_input <- opt[["sample_input"]]

#define a dataframe with centromere start and end
chromosome <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chr32", "chr33", "chr34", "chr35", "chr36", "chr37", "chr38", "chrX")

start.pos <- c(122678785,85426708,91889043,88276631,88915250,77573801,80974532,74330416,61074082,69331447,74389097,72498081,63241923,60966679,64190966,59632846,64289059,55844845,53741614,58134056,50858623,61439934,52294480,47698779,51628933,38964690,45876710,41182112,41845238,40214260,39895921,38810281,31377067,42124431,26524999,30810995,30902991,23914537,123869142)

end.pos <- c(122678785,85426708,91889043,88276631,88915250,77573801,80974532,74330416,61074082,69331447,74389097,72498081,63241923,60966679,64190966,59632846,64289059,55844845,53741614,58134056,50858623,61439934,52294480,47698779,51628933,38964690,45876710,41182112,41845238,40214260,39895921,38810281,31377067,42124431,26524999,30810995,30902991,23914537,123869142)

centromeres<-data.frame(chromosome=chromosome,start.pos=start.pos,end.pos=end.pos)

#function created using falcon to use the centromere position given
falcon.seg.seqz <- function(data.file, chromosome, centromeres){
    require(sequenza)
    seqz <- read.seqz(data.file, chr.name = chromosome)
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



data.file <- paste(input, "_small2.seqz.gz", sep="")
seqz.data <- read.seqz(data.file)
str(seqz.data, vec.len=2)

brk<-list()
chromosome.list <- c(1:38,"X")
for (i in chromosome.list){
    brk[[i]] <- falcon.seg.seqz(data.file, chromosome = i, centromeres = centromeres)
}
breaks <- do.call(rbind, brk)

gc.stats <- gc.sample.stats(data.file)
str(gc.stats)
options("scipen"=100, "digits"=4)
#this will force the program to read the regions as only straight values and not as exponential

gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]

chromosome.list <- c(1:38,"X")
sample <- sequenza.extract(sample_input, verbose=T,breaks=breaks , chromosome.list=chromosome.list)

#estimates cellularity and ploidy
CP<- sequenza.fit(sample, mc.core = 4, segment.filter= 5e6)

sequenza.results(sample, CP, out.dir=paste(input,"_falcon_seqz_centromeres",sep=""), sample.id = sample_name)

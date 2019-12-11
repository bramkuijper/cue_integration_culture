# select simulations that could be plotted

the.data <- read.table("../../data/summary_cue_int_posneg.csv"
        ,sep=";"
        ,header=T)

#the.data.s <- subset(the.data,qjuv==1.0 & qmat == 1.0 & p %in% c(0.1,0.9))
the.data.1 <- subset(the.data,qjuv==1.0 & qmat == 1.0 & p == 0.103)

stopifnot(nrow(the.data.1) > 0)

the.data.2 <- subset(the.data,qjuv==1.0 & qmat == 1.0 & p == 0.897)

stopifnot(nrow(the.data.2) > 0)

print(the.data.2)

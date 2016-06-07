argv <- commandArgs(TRUE)
infile <- argv[1]
outfile <- argv[2]

library(ggplot2)
library(reshape2)

pdf(outfile, width = 10, height = 6)
    fq <- read.table(infile, row.names = NULL, header = T, skip=3)
    t1 <- melt(fq[1:6], id = "pos")
    ggplot(t1, aes(pos, value, colour = variable)) + geom_line() + xlab("") + ylab("") + scale_colour_hue(name = "Base") + ylim(0,0.6)

    t2 <- fq[c(1, 7:49)]
    colnames(t2) <- c("pos", 0:42)
    t2 <- melt(t2, id = "pos")
    colnames(t2) <- c("pos", "Qual", "value")
    ggplot(t2, aes(x = pos, y = Qual)) + aes_string(fill = t2$value) + geom_tile() + scale_colour_gradient2() + scale_fill_gradient(low = "white", high = "purple", names("")) + scale_x_continuous(breaks = seq(0, 150, 10))
dev.off()

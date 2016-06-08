argv <- commandArgs(TRUE)
    infile1 <- argv[1]
    infile2 <- argv[2]
    outfile <- argv[3]

library(ggplot2)
library(reshape2)
library(gridExtra)

draw_distribution <- function(data, legend = TRUE) {
    data <- melt(data[1:6], id = "pos")
    p <- ggplot(data, aes(pos, value, colour = variable)) +
            geom_line() +
            labs(x = "", y = "") +
            scale_colour_hue(name = "Base") +
            ylim(0,0.6)
    if(legend == FALSE) {
        p <- p + theme(legend.position = "none") + 
                labs(title = "Reads 1")
    } else {
        p <- p + theme(axis.text.y = element_blank()) +
                labs(title = "Reads 2")
    }
    return(p)
}

draw_quality <- function(data, legend = TRUE) {
    data <- data[c(1, 7:49)]
    colnames(data) <- c("pos", 0:42)
    data <- melt(data, id = "pos")
    colnames(data) <- c("pos", "Qual", "value")

    p <- ggplot(data, aes(x = pos, y = Qual)) +
            aes_string(fill = data$value) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "purple", names("")) +
            scale_x_continuous(breaks = seq(0, 150, 10))
    if(legend == FALSE) {
        p <- p + theme(legend.position = "none") +
                labs(title = "Reads 1")
    } else {
        p <- p + theme(axis.text.y = element_blank()) +
                labs(x = "", y = "", title = "Reads 2")
    }
    return(p)
}

# read data
fq1 <- read.table(infile1, row.names = NULL, header = T, skip=3)
fq2 <- read.table(infile2, row.names = NULL, header = T, skip=3)
# draw plots
plot1 <- draw_distribution(fq1, FALSE)
plot2 <- draw_distribution(fq2)
plot3 <- draw_quality(fq1, FALSE)
plot4 <- draw_quality(fq2)
# arrangement && save
pdf(outfile, width = 12, height = 6)
    grid.arrange(plot1, plot2, ncol = 2)
    grid.arrange(plot3, plot4, ncol = 2)
dev.off()

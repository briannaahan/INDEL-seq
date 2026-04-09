install.packages("ggplot2")
install.packages("ggbeeswarm")
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
ID <- args[1]

path_distance <- paste("/Users/briannahan/Desktop/INDEL-seq-main/results/", ID, ".violin.csv", sep="")
distance <- read.csv(path_distance)
repeats <- subset(distance, isRepeat=="TRUE")
non_repeats <- subset(distance, isRepeat=="FALSE")

cpy1 <- data.frame(distance)
cpy2 <- data.frame(distance)
distance$group <- "All"
cpy1$group <- "Repeats"
cpy2$group <- "Non-Repeats"
distance <- rbind(distance,cpy1,cpy2)

repeats$group <- "Repeats"
non_repeats$group <- "Non-Repeats"


bind_rows(list(all = distance, repeats = repeats, non_repeats = non_repeats),
          .id = "Group") %>% 
  ggplot(aes(x=group, y=Distance, fill=Group)) + 
  geom_violin(position="identity", alpha=0.5,scale="count") +
  scale_y_continuous(trans='log10',labels = scales::comma) +
  scale_fill_manual(values=c("white", "red", "blue")) +
  xlab("") +
  ylab("Distance (bp)")


ggplot(data=distance, aes(x=NA, y=Distance)) +
  geom_violin() +
  geom_quasirandom(colour=distance$color, size = 2.5,alpha = 0.5,stroke=NA) +
  theme(legend.position = "none") +
  scale_y_continuous(trans='log10',labels = scales::comma)


#!/usr/bin/env Rscript
library(ggplot2)

#gene metaplot
print("Making combined gene metaplot")
df1 <- read.table("v2/results/gene_metaplot.tsv",header=T,sep="\t")
df2 <- read.table("v4/results/gene_metaplot.tsv",header=T,sep="\t")
df <- as.data.frame(cbind(df1,df2))
colnames(df) <- c("Bin","mCG1","mCHG1","mCHH1","Bin2","mCG2","mCHG2","mCHH2")
plot <- ggplot(df, aes(x=Bin)) + geom_line(aes(y=mCG1), color="dodgerblue4", size=0.8) +
              geom_line(aes(y=mCHG1), color="olivedrab", size=0.8) +
              geom_line(aes(y=mCHH1), color="hotpink4", size=0.8) +
              geom_line(aes(y=mCG2), color="dodgerblue2", size=0.8) +
              geom_line(aes(y=mCHG2), color="olivedrab2", size=0.8) +
              geom_line(aes(y=mCHH2), color="hotpink2", size=0.8) +
              theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
              axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
              legend.position="none", axis.line=element_line(color="black")) +
              ylab("Percent methylation") + xlab( "" ) +
              scale_y_continuous(limits=c(0,0.5), expand=c(0,0),
              breaks=c(0.15,0.30,0.45),labels=c("15%","30%","45%")) +
              geom_vline(xintercept=20, linetype="longdash", color="grey55") +
              geom_vline(xintercept=40, linetype="longdash", color="grey55") +
              scale_x_continuous(labels=c("-2000","TSS","TTS","+2000"), breaks=c(1, 20, 40, 60))
ggsave(filename="../figures_tables/combined_gene_metaplot.pdf", plot, height=4, width=4, useDingbats=F)

#gene metaplot
print("Making combined syntenic gene metaplot")
df1 <- read.table("v2/results/syntenic_gene_metaplot.tsv",header=T,sep="\t")
df2 <- read.table("v4/results/syntenic_gene_metaplot.tsv",header=T,sep="\t")
df <- as.data.frame(cbind(df1,df2))
colnames(df) <- c("Bin","mCG1","mCHG1","mCHH1","Bin2","mCG2","mCHG2","mCHH2")
plot <- ggplot(df, aes(x=Bin)) + geom_line(aes(y=mCG1), color="dodgerblue4", size=0.8) +
              geom_line(aes(y=mCHG1), color="olivedrab", size=0.8) +
              geom_line(aes(y=mCHH1), color="hotpink4", size=0.8) +
              geom_line(aes(y=mCG2), color="dodgerblue2", size=0.8) +
              geom_line(aes(y=mCHG2), color="olivedrab2", size=0.8) +
              geom_line(aes(y=mCHH2), color="hotpink2", size=0.8) +
              theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
              axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
              legend.position="none", axis.line=element_line(color="black")) +
              ylab("Percent methylation") + xlab( "" ) +
              scale_y_continuous(limits=c(0,0.5), expand=c(0,0),
              breaks=c(0.15,0.30,0.45),labels=c("15%","30%","45%")) +
              geom_vline(xintercept=20, linetype="longdash", color="grey55") +
              geom_vline(xintercept=40, linetype="longdash", color="grey55") +
              scale_x_continuous(labels=c("-2000","TSS","TTS","+2000"), breaks=c(1, 20, 40, 60))
ggsave(filename="../figures_tables/combined_syntenic_gene_metaplot.pdf", plot, height=4, width=4, useDingbats=F)

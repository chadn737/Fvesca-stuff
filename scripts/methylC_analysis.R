#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
#test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Sample name must be given", call.=FALSE)
}

#define plotting function for features
plot_features <- function(df){
  require(ggplot2)
  ggplot(df, aes(x=Bin)) + geom_line(aes(y=mCG), color="dodgerblue4", size=0.8) +
         geom_line(aes(y=mCHG), color="olivedrab", size=0.8) +
         geom_line(aes(y=mCHH), color="hotpink4", size=0.8) +
         theme(panel.background=element_blank(), panel.grid=element_blank(),
         axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
         axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
         legend.position="none", axis.line=element_line(color="black")) +
         ylab("Percent methylation") + xlab( "" ) +
         scale_y_continuous(limits=c(0,1), expand=c(0,0),
         breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
         geom_vline(xintercept=20, linetype="longdash", color="grey55") +
         geom_vline(xintercept=40, linetype="longdash", color="grey55")
}

#plot chr metaplots
if(file.exists("results/genome_windows_data.tsv")){
  print("Making genome distribution plots")
  df <- read.table("results/genome_windows_data.tsv",header=T,sep="\t")

  print("Making chromosome plots")
  library(stringi)
  df$chr <- stri_replace_last(df$window,'',regex="_[0-9]+.*?")
  df$window <- as.numeric(stri_replace(df$window,'',regex=".*_"))
  chr <- unique(df$chr)
  for(i in chr){
     tmp <- df[df$chr == i,]
     if (max(tmp$window) > 10){
       break_points <- c(2,4,6,8,10)*max(tmp$window)/10
       label_points <- c(2,4,6,8,10)*max(tmp$window)/20
       plot = ggplot(tmp, aes(x=window, group=1)) +
              geom_line(aes(y=mCG), color="dodgerblue4", size=0.8) +
              geom_line(aes(y=mCHG), color="olivedrab", size=0.8) +
              geom_line(aes(y=mCHH), color="hotpink4", size=0.8) +
              theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
              axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
              legend.position="none", axis.line=element_line(color="black")) +
              ylab("Percent methylation") + xlab( "Kbps" ) +
              scale_y_continuous(limits=c(0,1), expand=c(0,0),
              breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
              scale_x_continuous(expand=c(0,0), breaks=break_points, labels=label_points)
       filename=paste("figures_tables/Fvesca_", args[1], "_Chr", i,"_metaplot.pdf", sep="")
       ggsave(filename=filename, plot, height=4, width=8, useDingbats=F)
       rm(tmp,plot)
     }
  }
  rm(df)
}

#plot gene metaplots
if(file.exists("results/gene_metaplot.tsv")){
  print("Making gene metaplot")
  df <- read.table("results/gene_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/Fvesca_", args[1], "_gene_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}

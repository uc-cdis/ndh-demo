#install.packages("IRkernel")
#install.packages('dplyr')
#library(dplyr)
#install.packages("htmlwidgets")
#install.packages("ggplot2")
#install.packages('pheatmap')
#install.packages("devtools")
library(devtools)
#install_github("hrbrmstr/waffle")
library(ggplot2)
library(waffle)
library(htmlwidgets)
library(pheatmap)
library(grid)
library(gridExtra)

heat_map <- function(){
  clean_df <- read.csv('out.summary.csv',header = TRUE, stringsAsFactors = FALSE)
  clean_df[clean_df == 'yes'] <- 1
  clean_df[clean_df == 'no'] <- 0
  rownames(clean_df) <- clean_df$name
  plot_df <- data.matrix(clean_df[,2:57])
  hmcols<- colorRampPalette(c("pink","lawngreen"))(2)

  p =pheatmap(plot_df, cluster_row=T,cluster_cols =F,color=hmcols,
     fontsize=10, fontsize_col=5,scale="none",legend_breaks=c(0,0.25,0.5,0.75,1),
        legend_labels=c('','Not match with cluster','','Match with cluster',''),
         show_colnames = T,show_rownames=F)
  p
}

scatter_plot <- function(){
  d=data.frame(Sensitivity=c(1,0.6803,0.494,0.975,0.947,0.724,0.988,0.95),
   Specificity=c(0.831,0.508,0.975,0.248,0.924,0.859,0.953,0.61),
   Method=c('Ariba','Ariba','Ariba','Ariba','Mykrobe','Mykrobe','Mykrobe','Mykrobe'),
   Drug=c('Rifampicin','Ethambutol','Isoniazid','Pyrazinamide','Rifampicin','Ethambutol','Isoniazid','Pyrazinamide'),
   info=c('Rifampicin_Ariba','Ethambutol_Ariba','Isoniazid_Ariba','Pyrazinamide_Ariba','Rifampicin_Mykrobe','Ethambutol_Mykrobe','Isoniazid_Mykrobe','Pyrazinamide_Mykrobe'))
  ggplot()  + xlim(0, 1) + ylim(0, 1) +
   geom_point(data=d, mapping=aes(x=Specificity, y=Sensitivity, color=Drug, shape=Method),size=4) +
   ggtitle("Prediction Tools Comparison") + scale_color_brewer(palette="Set1")+
    theme(plot.title = element_text(hjust = 0.5))+
   geom_text(aes(label = d$info),size = 2,x=d$Specificity-0.088, y=d$Sensitivity)+ theme_minimal() +
   theme(legend.position = c(0.2, 0.4))
}

grid_plot <- function(){
  #A-RIFA
  p1 <- waffle(c(Match = 600, Not_Match = 56,Not_Available=0), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
  #A-ETH
  p2 <- waffle(c(Match = 413, Not_Match = 243,Not_Available=0), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
#A-ISO
  p3 <- waffle(c(Match = 519, Not_Match = 137,Not_Available=0), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
#A-PYR
  p4 <- waffle(c(Match = 539, Not_Match = 113, Not_Available=4), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
#M-RIFA
  p5 <- waffle(c(Match = 613, Not_Match = 43,Not_Available=0), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
#M-ETH
  p6 <- waffle(c(Match = 502, Not_Match = 154,Not_Available=0), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
#M-ISO
  p7 <- waffle(c(Match = 634, Not_Match = 22,Not_Available=0), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))
#M-PYR
  p8 <- waffle(c(Match = 578, Not_Match = 74,Not_Available=4), rows = 16, keep = FALSE,size = .2,colors = c("lawngreen","pink",'grey'),
  legend_pos = c(16,16))

  col1 <- arrangeGrob(p2, p3, p4,p1, top="Ariba",ncol=1)
  col2 <- arrangeGrob(p6,p7,p8,p5, top='Mykrobe',ncol=1)
  col3 <- textGrob(x = 0.5, y = 0.48, paste("Ethambutol\n\n\n\n\n", "Isoniazid\n\n\n\n",
        "\nPyrazinamide\n\n\n\n\n", "Rifampicin"),gp=gpar(fontsize=9,font=8))
  legd <- legendGrob(c("Match with phenotype", "Not match with phenotype", "Not Available"), pch=22,
   nrow=1, byrow=TRUE, gp=gpar(col = c("lawngreen","pink",'grey'), fill = c("lawngreen","pink",'grey')))
  grid.arrange(col1,col3,col2, legd, widths=c(1,0.2,1), heights=c(4,0.5),top = "Accuracy Comparison",
             layout_matrix=rbind(c(1,2,3), c(4,4,4)))
}


sorted_grid <- function() {
    t=read.csv('match_plot_sort_mycrobe.txt',sep='\t')
 colors <- colorRampPalette(c("lawngreen", "grey", "pink"))(3)
 ggplot(t[1:5248,], aes(column, row, fill = match_not)) + 
  geom_tile(colour = "white") + 
  facet_grid(drug~method)  +theme_minimal()+
  labs(x="", y="") + ggtitle("Side by Side Comparison")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values=colors)+  theme(strip.text.x = element_text(size = c(10)))+
  theme(strip.text.y = element_text(size = c(7)))
}

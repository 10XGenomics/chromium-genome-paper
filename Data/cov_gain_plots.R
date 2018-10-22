library(TenXploreR)
library(ggplot2)
library(reshape)
library(cowplot)
library(dplyr)
library(wesanderson)

#Read in the results from the cov_gain_calculator.py script
result <- read.csv("cov_gain_results.csv")
result$sequencer <- as.factor(result$sequencer)
result$Var1_length_sweep <- as.factor(result$id1_length_sweep)
result <- result[,c("cell_line","id1_length_sweep","Cov_in_10x","Cov_in_TruSeq","Delta")]
result <- melt(result)

#Functions to make main boxplots
boxplot <- function(df, color_by, x_axis_labels,legend_title,color){
    ggplot(df, aes(x=variable,y=value/1000000)) +
    geom_boxplot(outlier.color = "white") +
    geom_point(aes_string(colour=color_by),position=position_jitter(width=0.2, height=0),size=2) +
    scale_color_manual(name=legend_title,values=color) +
    theme(plot.margin = unit(c(20,20,20,20),"pt")) +
    scale_x_discrete(labels=x_axis_labels) + ylab("Mb") +xlab("") +
    theme(axis.text.x=black.10.45.text, axis.text.y=element_text(color="black",size=14)) +
    theme(legend.text=element_text(color="black",size=14))
}

boxplot_noleg <- function(df, color_by, x_axis_labels,legend_title,color){
  ggplot(df, aes(x=variable,y=value/1000000)) +
    geom_boxplot(outlier.color = "white") +
    geom_point(aes_string(colour=color_by),position=position_jitter(width=0.2, height=0),size=2) +
    scale_color_manual(name=legend_title,values=color) +
    theme(plot.margin = unit(c(20,8,8,20),"pt")) +
    scale_x_discrete(labels=x_axis_labels) + ylab("Mb") +xlab("") +
    theme(axis.text.x=black.10.45.text, axis.text.y=element_text(color="black",size=14)) +
    theme(legend.position="none",legend.text=element_text(color="black",size=14))
}

x_axis_labels = c("Covered in 10x \n not in TruSeq","Covered in TruSeq \n not in 10x","Net coverage \n gain in 10x")
black.10.45.text <- element_text(color = "black", size = 14, angle= 45, hjust = 1)

plot1 <- boxplot(filter(result,id1_length_sweep=="No"),"cell_line",x_axis_labels,"Cell lines",wes_palette("GrandBudapest1",3,type="discrete"))
save_plot("Fig2_Coverage_Gain.png",plot1,ncol=1,nrow=1,base_aspect_ratio = 1.3)

plot2 <- boxplot_noleg(filter(result,id1_length_sweep!="No"),"id1_length_sweep",x_axis_labels,"Molecule length",wes_palette("Zissou1",5))
inset1_df <- filter(result,id1_length_sweep!="No" & variable=="Cov_in_10x")
inset2_df <- filter(result,id1_length_sweep!="No" & variable=="Cov_in_TruSeq")
inset1 <-boxplot_noleg(inset1_df,"id1_length_sweep",x_axis_labels[1],"Molecule length",wes_palette("Zissou1",5))
inset2 <-boxplot(inset2_df,"id1_length_sweep",x_axis_labels[2],"Molecule length",wes_palette("Zissou1",5))

ls.plot <- plot_grid(plot2,inset1,inset2,ncol=3,align="h",labels=c("A","B","C"),rel_widths=c(1,1,1.5))
save_plot("SuppFig3_Coverage_Gain_LengthSweep.png",ls.plot,ncol=3,base_aspect_ratio=0.8)




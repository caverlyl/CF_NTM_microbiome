#Code adapted from Ridenhour, Benjamin J., et al. "Modeling time-series data from microbial communities." The ISME journal 11.11 (2017): 2526-2537.

source("data_prep_dis_all.R")
library(corrplot)
library(fields)
library(reshape2)
library(igraph)
library(dplyr)
library(funModeling)
results <- read.csv("PoissonARIMA-Disease_dis_all.csv")

results <- results[results$i !=1, ]

results <- results[results$dev.ratio > 0.02,] #drop models will low deviance ratios
meanReads <- colMeans(X)
meanReads[1:23] <- 1 #subject effects need to have a 1
meanReads <- c(1,meanReads) #add a 1 for the intercept
results$x.scaled <- results$x * meanReads[results$i] #create scaled coefficients

results[,1] <- append("Intercept",colnames(X))[results[,1]]
results[,2] <- colnames(Y)[results[,2]]

df <- results[,c("i","j","x.scaled")]

df[,c("i")] <- sub("as.factor\\(oxy.full\\$Subject\\)","Patient_",df[,c("i")])

mat <- xtabs(df[, 3] ~ df[, 2] + df[, 1])
mat[mat == 0] <- NA

#add taxa names
taxmap <- read.csv("taxa_key.csv",header=FALSE,row.names=1)
old_mat <- mat
rownames(mat)<-taxmap[rownames(mat),]
colnames(mat)<-taxmap[colnames(mat),]

#Heatmap
pdf("Disease_all_coef_matrix.pdf")
corrplot(as.matrix(mat[,c(1:16)]),addgrid.col = "white",tl.col = "black",  na.label = "square", na.label.col = "gray",is.corr=FALSE,tl.cex=0.6,number.cex=1,method="color")
dev.off()

#Network
input <- as.matrix(mat[,c(1:16)])
my_cor_df <- melt(t(input))
colnames(my_cor_df) <- c("Var1","Var2","value")
my_cor_df <- filter(my_cor_df, value != "NA")
my_cor_df$Var1 <- as.character(my_cor_df$Var1)
my_cor_df$Var2 <- as.character(my_cor_df$Var2)
my_cor_df <- filter(my_cor_df, Var1 != Var2)

net <- graph.data.frame(my_cor_df, directed = TRUE)
net1 <- simplify(net,remove.loops = FALSE)
ceb <- cluster_optimal(net1)

c_scale <- colorRamp(c('red',"lightgray",'blue'))
weights <- append(my_cor_df[,3],c(-1*max(my_cor_df[,3]),-1*min(my_cor_df[,3])))
E(net)$color = apply(c_scale(range01(weights))[1:28,], 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
ordered_colors = apply(c_scale(range01(seq(from=min(weights),to=max(weights),by=0.001))), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )

tester = membership(ceb)

sums <- as.data.frame(meanReads)
sums[,2] <- taxmap[row.names(sums),]
names = rownames(as.data.frame(net[1]))
size <- as.data.frame(names)
row.names(size) <- size[,1]
for (i in 1:length(size$names)){

    size$s[i] = sum(sums[sums[,2] == names[i],][1],na.rm=TRUE)

}
sizes <- as.data.frame(size$s)
rownames(sizes)<- rownames(size)




pdf(file="NTM_disease_all_network.pdf")
plot(net, edge.arrow.size=0.5, edge.width = range01(abs(my_cor_df[,3])), vertex.color=tester, vertex.frame.color=NA, vertex.label.cex=0.7,vertex.label.dist=1.5,vertex.size=20*sizes$s^(1/3),vertex.label.degree=pi/2,vertex.label.color="black")
image.plot( legend.only=TRUE, axis.args=list(cex.axis=0.6), zlim= c(min(weights),max(weights)),col=ordered_colors,smallplot= c(.05,.065, .2, .4))
dev.off()
#Code adapted from Ridenhour, Benjamin J., et al. "Modeling time-series data from microbial communities." The ISME journal 11.11 (2017): 2526-2537.

mb <- read.csv("interpolated_counts_3_5_with_metadata_fev1_dis_yes.csv")

#the way the data are being processed for the AR(1)
#utilizes the fact that row i and row i+1 represent a samples
#from t and t+1 respectively

#The next 4 lines interject an extra row between
#time series from different individuals. This helps with formatting
#the data later; specifically these rows will be used to eliminate
#data where we did not have samples at time t and t+1 for an
#individual.
extraStep <- c(unique(mb$Time.Point),1000)
extraSet2<-expand.grid(extraStep,unique(mb$Animal))
names(extraSet2)<-c("Time.Point","Animal")
oxy.full <-merge(extraSet2,mb,all.x=T) #merge in the fake timepoints, note these all show up as rows of NAs

#sort the data appropriately for forming the time series
oxy.full <- oxy.full[do.call("order",oxy.full[,c("Animal","Time.Point")]),]

#create a model matrix for the individual factor (Animal ID)
animal.matrix<-as.matrix(model.matrix(~as.factor(oxy.full$Animal)-1))

#drop some of the metadata that will not be used for the analysis
#oxy.full <- oxy.full[, ! names(oxy.full) %in% c("Time.Point","Animal","Sample","Oxalate.Concentration")]

oxy.full <- oxy.full[, ! names(oxy.full) %in% c("Time.Point","Animal")]

#Because of the large number of OTUs, many of which are largely 0s,
#we need to filter the data some. The top 10% will be kept (data frame top10).
#otuReads <- colSums(oxy.full[,! names(oxy.full) %in% c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads")], na.rm = T) #sum the reads for each OTU
#keepOTU <- names(otuReads)[otuReads > quantile(otuReads,0.9)] #get the names of the OTUs in the top 10%
#top10 <- oxy.full[,union(keepOTU,c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads"))]

top10 <- oxy.full

#a quick histogram of the total read counts for the OTUs in 
#hist(log(colSums(top10[,! names(top10) %in% c("NTM_disease","Persistent_infection", "Total.Reads")],na.rm=T)), xlab = "Log(OTU Reads)",main = NULL)

#create a couple of vectors indicating which columns should be
#in the X and Y matrices
xNames <- setdiff(names(top10), c("fev1","NTM_disease","Persistent_infection", "Total.Reads")) #gives strictly OTU names (i.e. denovoXXXXXX)
yNames <- setdiff(names(top10), c("NTM_disease","Persistent_infection")) #gives all OTUs in top10 + oxalate excreted and degraded

#create the X matrix (time t) and Y matrix (time t+1)
X <- top10[1:(nrow(top10)-1),]
Y <- top10[2:nrow(top10), ]

X[,xNames] <- X[,xNames] / X[,"Total.Reads"] #make OTU data in X relative abundances
#X <- X[,!names(X) %in% c("NTM_disease","Persistent_infection")]  #drop excreted and degraded
X <- cbind(animal.matrix[1:(nrow(animal.matrix)-1),], X) #add in the model matrix for subject offsets

#the next command uses NAs to identify which rows of X and Y should
#be dropped (because of missing data); any rows where oxalated consumed
#is NA in either X or Y indicates a bad data point
useMe<- ! (is.na(Y$NTM_disease)|is.na(X$NTM_disease))
X <- X[useMe,]
Y <- Y[useMe,]

#because of the study design, oxalate consumption actually needs
#to be shifted (i.e. have a different lag)
#X$Oxalate.Consumed <- Y$Oxalate.Consumed 

X <- X[,!names(X) %in% c("Total.Reads","Persistent_infection","NTM_disease","fev1")]
Y <- Y[,!names(Y) %in% c("NTM_disease","Persistent_infection","Total.Reads","fev1")]

#write to the global environment
goodOTUs  <- colMeans(X) > 0.01
goodOTUs[1:8] <- TRUE
X <- X[,goodOTUs]
Y <- Y[,goodOTUs[9:length(goodOTUs)]]
X <<- X
Y <<- Y
#oxy.taxa <<- oxy.taxa
 #end with()

#cleanup
#rm(e)

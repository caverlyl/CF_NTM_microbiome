
#Code adapted from Ridenhour, Benjamin J., et al. "Modeling time-series data from microbial communities." The ISME journal 11.11 (2017): 2526-2537.

mb <- read.csv("interpolated_counts_3_5_with_metadata_fev1_dis_all.csv")

#the way the data are being processed for the AR(1)
#utilizes the fact that row i and row i+1 represent a samples
#from t and t+1 respectively

#The next 4 lines interject an extra row between
#time series from different individuals. This helps with formatting
#the data later; specifically these rows will be used to eliminate
#data where we did not have samples at time t and t+1 for an
#individual.
extraStep <- c(unique(mb$Time.Point),1000)
extraSet2<-expand.grid(extraStep,unique(mb$Subject))
names(extraSet2)<-c("Time.Point","Subject")
oxy.full <-merge(extraSet2,mb,all.x=T) #merge in the fake timepoints, note these all show up as rows of NAs

#sort the data appropriately for forming the time series
oxy.full <- oxy.full[do.call("order",oxy.full[,c("Subject","Time.Point")]),]

#create a model matrix for the individual factor (Subject ID)
Subject.matrix<-as.matrix(model.matrix(~as.factor(oxy.full$Subject)-1))

#drop some of the metadata that will not be used for the analysis
oxy.full <- oxy.full[, ! names(oxy.full) %in% c("Time.Point","Subject")]


#create a couple of vectors indicating which columns should be
#in the X and Y matrices
xNames <- setdiff(names(oxy.full), c("fev1","NTM_disease","Persistent_infection", "Total.Reads")) #gives strictly OTU names (i.e. denovoXXXXXX)
yNames <- setdiff(names(oxy.full), c("NTM_disease","Persistent_infection")) 

#create the X matrix (time t) and Y matrix (time t+1)
X <- oxy.full[1:(nrow(oxy.full)-1),]
Y <- oxy.full[2:nrow(oxy.full), ]

X[,xNames] <- X[,xNames] / X[,"Total.Reads"] #make OTU data in X relative abundances
X <- cbind(Subject.matrix[1:(nrow(Subject.matrix)-1),], X) #add in the model matrix for subject offsets

#the next command uses NAs to identify which rows of X and Y should
#be dropped (because of missing data)
useMe<- ! (is.na(Y$NTM_disease)|is.na(X$NTM_disease))
X <- X[useMe,]
Y <- Y[useMe,]

X <- X[,!names(X) %in% c("Total.Reads","Persistent_infection","NTM_disease","fev1")]
Y <- Y[,!names(Y) %in% c("NTM_disease","Persistent_infection","Total.Reads","fev1")]

#filter OTUs by mean relative abundance greather than 0.01 
goodOTUs  <- colMeans(X) > 0.01
goodOTUs[1:23] <- TRUE
X <- X[,goodOTUs]
Y <- Y[,goodOTUs[24:length(goodOTUs)]]

#write to the global environment
X <<- X
Y <<- Y


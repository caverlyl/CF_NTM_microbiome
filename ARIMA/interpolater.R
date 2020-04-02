
#Read counts scaled to ~1660 for each sample prior to interpolation
df <- read.csv(file="counts_3_5_scaled_fev1.csv",row.names=1)


all_patients=c()
for(i in unique(df[,c("Patient_ID")])){
    results = c()
    labels = c()
    current = df[df[,"Patient_ID"] ==i,]
    min_age = min(current[,c("relative_age")])
    max_age = max(current[,c("relative_age")])
    c=ceiling(2*min_age)/2
    f=floor(2*max_age)/2
    ticks = seq(from=c,to=f,by=0.5)
    x <- current[,c("relative_age")]
    for(j in colnames(df)[5:49]){
        y <- current[,j]
        interp = approx(x,y,ticks)
        results <- cbind(results,interp$y)
        labels <- cbind(labels,j)

    }
    results <-cbind(results,ticks)
    results <-cbind(results,rep(i,length(ticks)))
    results <-cbind(results,rep(df[df[,"Patient_ID"] ==i,]$NTM_disease[1],length(ticks)))
    results <-cbind(results,rep(df[df[,"Patient_ID"] ==i,]$Persistent_infection[1],length(ticks)))
    labels <- cbind(labels,"relative_age","Patient_ID","NTM_disease","Persistent_infection")
    colnames(results)<-labels
    all_patients <- rbind(all_patients,results)

}
rounded <- round(all_patients[,1:44])
final <- cbind(rounded,all_patients[,45:49])

write.csv(final,file="interpolated_counts_3_5_with_metadata_fev1.csv")
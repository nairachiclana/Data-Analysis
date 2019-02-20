options(warn=-1)

load(file="savedData/filtered_df.RData")

last.time<-max(na.omit(as.numeric(filtered.df$time)))
replicated.df<-data.frame()

total=nrow(filtered.df)
pb <- txtProgressBar(min = 0, max = total, style = 3)

for (patient in 1:nrow(filtered.df)) {
  setTxtProgressBar(pb, patient)
  time<-filtered.df$time[patient]
  status<-filtered.df$status[patient]
  patient.row<-filtered.df[patient,]
  
  #time before patient time->status=0
  for(t in seq(0,time,by=20)) {
    patient.row$time<-t
    patient.row$status<-0
    replicated.df<-rbind(replicated.df, patient.row)
  }
  #if deceased->time after deceased->status=1
  if(status==1) {
    for(t in seq(time,last.time,by=20)) {
      patient.row$time<-t
      patient.row$status<-1
      replicated.df<-rbind(replicated.df, patient.row)
    }
  }
}

names(replicated.df)<-colnames(filtered.df)

cat("\n", "Las dimensiones del dataframe replicado son:", dim(replicated.df), "\n")

save(replicated.df, file="savedData/replicated_df.RData")

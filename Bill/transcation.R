#####constant definition#######
# s=0
# w=0
# ya=0
# yu=0
# z=0
# 
# names(s)="suyang"
# names(w)="wang"
# names(ya)="yang"
# names(yu)="yuliang"
# names(z)="zhang"
# 
# family=c(s,w,ya,yu,z)

#####main##################

account <- function(amount , payer, except=c()){
  
  transaction <- function(amount , payer, except=c()){
    i=1
    Nexc=length(except)
    Locpayer=which(names(family)==names(payer))
    participant=family
    while (i<=Nexc){
      temp=which(names(participant)==names(except[i]))
      participant=participant[-temp]
      i=i+1
    }
    Npart=length(participant)
    family[Locpayer]=family[Locpayer]+amount/Npart*(Npart-1)
    temp=which(names(participant)==names(payer))  #loction of payer in participant
    vicepart=participant[-temp]
    for(i in 1:Npart-1){
      tempLoc=which(names(family)==names(vicepart[i]))
      family[tempLoc]=family[tempLoc]-amount/Npart
    }
    details=c(paste(names(participant),collapse=","),amount/Npart,names(payer),amount/Npart*(Npart-1),amount,format(Sys.Date(), "%Y.%m.%d"))
    names(details)=c("participant","share","payer","PayerGain","total","date")
    result <- list(family, data.frame(details))
    names(result)<-c("family","details")
    return(result)
  }

details=read.csv("details.csv",header=TRUE)[-1]
family=read.csv("family.csv",header=TRUE)[-1]

this<-transaction(amount,payer,except)

temp=rbind(details,t(this$details))
write.csv(temp,file="details.csv")
write.csv(t(this$family),file="family.csv")

}
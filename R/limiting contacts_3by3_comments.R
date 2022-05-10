#this code creates a 3-by-3 "ready reckoner"
#accompanying the paper Brooks-Pollock et al. 2021 
#Mapping social distancing measures to the reproduction number for COVID-19
#https://doi.org/10.1098/rstb.2020.0276
#Code developed by Ellen Brooks-Pollock 2020/2021

#load necessary libraries
library(tidyverse)
library(boot)

# DATA SET UP ---------------------------------

#1. load Social Contact Survey data
#2. create additional variables

#load processed Social Contact Survey data 
#raw data can be downloaded from http://wrap.warwick.ac.uk/54273/
#Citations for data 
#Danon et al. 2012 doi:10.1098/rsif.2012.0357
#Danon et al. 2013  doi:10.1098/rspb.2013.1037
#the data contain one line per contact reported
#PID is the person ID of the repondent reporting the contact
load(file="data/contactsperson_withkids6.Rdata")


#determine probability individual has symptoms if infected
#I assume the probability of symptoms increases linearly with age from
#25% in the youngest age group to 75% in the oldest age groups
a1=1:100
hassymp = rep(0.25,length(a1))
hassymp[80:100] = 0.75
hassymp[19:80] = 0.25 + 0.5 *(a1[19:80]-a1[19])/(a1[80]-a1[19])
contactsperson$hassymp = 0*contactsperson$age
for(a1 in 1:length(contactsperson$age))
{
  if(!is.na(contactsperson$age[a1]))
  {
    if(contactsperson$age[a1]>=0)
    {
      contactsperson$hassymp[a1] = hassymp[contactsperson$age[a1]+1]
    }
  }
}

#define variables needed later
contactsperson$comply = 0*contactsperson$age
contactsperson$contacttrace = rep(FALSE,length(contactsperson$age))



# FUNCTIONS ---------------------------------

#function for population-level reproduction number
myweightedmean=function(DD,ix)
{
  return(sum(DD[ix,1]^2*DD[ix,2])/sum(DD[ix,2]))
}

#degreeperson function takes the list of contacts and groups by PID (person ID) to calculate degree and weighted degree
degreeperson = function(contactsperson,inf_child,CTF1,covidsec)
{
  contactsperson %>%
    filter(dcat>=0)%>%
    mutate(trans=ifelse(age<=11 | Home==1,1,covidsec)) %>%
    mutate(under18 = ifelse(age<=18,1,0)) %>%
    mutate(infect = ifelse(age<=11,inf_child,1*trans)) %>%
    mutate(Num1 = ifelse(contacttrace==TRUE,((1-CTF1)*Numbers)*(1-comply)*infect,Numbers*(1-comply)*infect)) %>%
    mutate(homenum = Num1*Home, worknum=Num1*WorkSchool, 
           othernum=Num1*OtherLeisure+Num1*Travel)%>%
    group_by(PID)%>%
    summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
              totdcat=sum(dcat),
              wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),ageweight=mean(ageweight)) -> degred
  
  return(degred)
  
}

#degreeperson0 function is similar to degreeperson but is for the first call only
degreeperson0 = function(contactsperson,inf_child,CTF1)
{
  contactsperson %>%
    filter(dcat>=0)%>%
    mutate(under18 = ifelse(age<=18,1,0)) %>%
    mutate(infect = ifelse(age<=11,inf_child,1)) %>%
    mutate(Num1 = ifelse(contacttrace==TRUE,((1-CTF1)*Numbers)*(1-comply)*infect,Numbers*(1-comply)*infect)) %>%
    mutate(homenum = Num1*Home, worknum=Num1*WorkSchool, 
           othernum=Num1*OtherLeisure+Num1*Travel)%>%
    group_by(PID)%>%
    summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
              totdcat=sum(dcat),
              wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),ageweight=mean(ageweight)) -> degred
  
  return(degred)
  
}


# CALIBRATION to baseline  ---------------------------------

# calibrate to no social distancing

#relative infectiousness of children under 11
infchild=1

#call degreeperson0 once to calculate degree per person
degall=degreeperson0(contactsperson,infchild,0.0)



#starting average R0
meanR0=2.8

#number of bootstrap replicates (reduce to 100 increase running speed during tests, 1000 for final figures)
NR=100
#NR=1000

#bootstrapped confidence interval required 
myinterval=0.95

#calculate the multiplier to give a popoulation-level reproduction number of meanR0
mult1=meanR0/(sum((degall$wcon^2)*degall$ageweight)/sum(degall$ageweight))

#test baseline values:
#rind=cbind(degall$wcon,degall$ageweight)
#myboot=boot(rind,statistic=myweightedmean,R=500)
#mult1*median(myboot$t)
#bci=boot.ci(myboot,type="perc",conf=myinterval)
#mult1*bci$percent[1,4:5]

rind=mult1*degall$wcon
hist(rind)

# MAIN RECKONER CALCULATIONS ---------------------------------


#set up colours for figures
mycols2=c(rgb(0,0.5,1,0.5),rgb(0,0.8,0.8,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.3,0.3,0.3,0.6),rgb(0.6,0.0,0.6,0.6))
par(mar=rep(1.7,4),mfrow=c(3,3),oma=c(9,4,2,4),xpd=FALSE)

#COMP is values of compliance between 
#0 (none compliant, i.e. contact takes place) and 
#1 (compliant with social distancing and contact does not take place)
COMP=seq(0,1,0.1)

#here we create a 3-by-3 figure for 3 levels of COVID-security and 3 levels of contact tracing
for(covidsec in c(1,0.75,0.5))
{
  x=1
  Rt=matrix(0,nrow=length(COMP),ncol=3)

  Rt5611_noCT=matrix(0,nrow=length(COMP),ncol=3)
  Rtprimary_noCT=matrix(0,nrow=length(COMP),ncol=3)
  Rtall_noCT=matrix(0,nrow=length(COMP),ncol=3)
  
  Rtclosed_20CT=matrix(0,nrow=length(COMP),ncol=3)
  Rt5611_20CT=matrix(0,nrow=length(COMP),ncol=3)
  Rtprimary_20CT=matrix(0,nrow=length(COMP),ncol=3)
  Rtall_20CT=matrix(0,nrow=length(COMP),ncol=3)
  
  Rtclosed_80CT=matrix(0,nrow=length(COMP),ncol=3)
  Rt5611_80CT=matrix(0,nrow=length(COMP),ncol=3)
  Rtprimary_80CT=matrix(0,nrow=length(COMP),ncol=3)
  Rtall_80CT=matrix(0,nrow=length(COMP),ncol=3)
  
  numreps=1
  x=1
  MAXCON=15
  for(cc in COMP)
  {
    
    #get relevant indices of different contact types
    ixO18=which(contactsperson$age>18) #all contacts of persons over 18
    ixU18=which(contactsperson$age<=18 & contactsperson$WorkSchool==1) #school contacts of persons <= 18
    ixU18_O=which(contactsperson$age<=18 & contactsperson$WorkSchool==0) #non-school contacts of persons <= 18
    ixO18W=which(contactsperson$age>18 & contactsperson$WorkSchool==1) #work contacts of persons > 18
    ixO18L=which(contactsperson$age>18 & contactsperson$WorkSchool==0) #non-work contacts of persons > 18
    ixH = which(contactsperson$Home==1) # home contacts
    
    #set up arrays for individual reproduction numbers
    
    #no contact tracing
    rind=c() #schools closed
    rind_5611=c()  #transition primary years back at school
    rind_primary=c() #primary schools open
    rind_allschool=c() #primary and secondary schools open
    
    #20% efficacy contact tracing 
    rind_closed_20CT=c() #schools closed
    rind_5611_20CT=c() #transition primary years back at school
    rind_primary_20CT=c() #primary schools open
    rind_allschool_20CT=c() #primary and secondary schools open
    
    #60% efficacy contact tracing 
    rind_closed_80CT=c() #schools closed
    rind_5611_80CT=c() #transition primary years back at school
    rind_primary_80CT=c() #primary schools open
    rind_allschool_80CT=c() #primary and secondary schools open
    
    #repeat random draws numreps times. Helps to reduce noise if number of samples is small 
    #numreps=5 is sufficient
    for(r in 1:numreps) 
    {
      #for each iteration, draw a random number to see if contact tracing took place
      contactsperson$contacttrace = 0*contactsperson$dcat
      myrands = runif(length(contactsperson$hassymp),min=0,max=1)
      contactsperson$contacttrace = myrands < contactsperson$hassymp
      #contactsperson$masstest = rep(FALSE,length(contactsperson$age))
      
      ####now: schools closed
      #for each reported contact, sample which contacts are not active (1=not active)
      contactsperson$comply = 0*contactsperson$dcat
      contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
      #98% of school contacts do not take place: 
      contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
      contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
      contactsperson$comply[ixH]=0 #home contacts
      
      ########NO CONTACT TRACING
      #calculate degree per person
      degred=degreeperson(contactsperson,infchild,CTF1=0.0,covidsec)
      #and the individual-level reproduction numbers
      rind=rbind(rind,cbind(degred$wcon,degred$ageweight))  
      
      ########60% CONTACT TRACING
      Lix=which(contactsperson$comply==0)
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.6)
      #CTF1=0.6
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_closed_80CT=rbind(rind_closed_80CT,cbind(degred$wcon,degred$ageweight))
      
      ########20% CONTACT TRACING
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.2)
      #CTF1=0.2
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_closed_20CT=rbind(rind_closed_20CT,cbind(degred$wcon,degred$ageweight))
      
    
      #ages 5, 6, 11 back
      contactsperson$comply = 0*contactsperson$dcat
      contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
      
      #98% of school contacts do not take place
      contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
      contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
      contactsperson$comply[ixH]=0 #home contacts
      
      #re-instate transition years:
      ixtrans=which((contactsperson$age==11 | contactsperson$age==5 | contactsperson$age==6 | contactsperson$age==7| contactsperson$age==8) & contactsperson$WorkSchool==1)
      contactsperson$comply[ixtrans]=0

      ########NO CONTACT TRACING, 20% CONTACT TRACING, 60% CONTACT TRACING
      degred=degreeperson(contactsperson,infchild,0.0,covidsec)
      rind_5611=rbind(rind_5611,cbind(degred$wcon,degred$ageweight))
      
      Lix=which(contactsperson$comply==0)
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.6)
      #CTF=0.6
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_5611_80CT=rbind(rind_5611_80CT,cbind(degred$wcon,degred$ageweight))
      
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.2)
      #CTF=0.2
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_5611_20CT=rbind(rind_5611_20CT,cbind(degred$wcon,degred$ageweight))
      
      
      #primary schools back
      contactsperson$comply = 0*contactsperson$dcat
      contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
      contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
      contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
      contactsperson$comply[ixH]=0 #home contacts
      
      #reinstate primary school aged school contacts:
      ixtrans=which((contactsperson$age<=11 & contactsperson$age>=5) & contactsperson$WorkSchool==1)
      contactsperson$comply[ixtrans]=0
      
      ########NO CONTACT TRACING, 20% CONTACT TRACING, 60% CONTACT TRACING
      degred=degreeperson(contactsperson,infchild,0.0,covidsec)
      rind_primary=rbind(rind_primary,cbind(degred$wcon,degred$ageweight))
      
      Lix=which(contactsperson$comply==0)
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.6)
      #CTF1=0.6
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_primary_80CT=rbind(rind_primary_80CT,cbind(degred$wcon,degred$ageweight))
      
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.2)
      #CTF1=0.2
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_primary_20CT=rbind(rind_primary_20CT,cbind(degred$wcon,degred$ageweight))
      
      #all schools back
      contactsperson$comply = 0*contactsperson$dcat
      contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
      contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
      contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
      contactsperson$comply[ixH]=0 #home contacts
      
      #reinstate all school contacts:
      ixtrans=which((contactsperson$age<=18 & contactsperson$age>=5) & contactsperson$WorkSchool==1)
      contactsperson$comply[ixtrans]=0
      
      ########NO CONTACT TRACING, 20% CONTACT TRACING, 60% CONTACT TRACING
      degred=degreeperson(contactsperson,infchild,0.0,covidsec)
      rind_allschool=rbind(rind_allschool,cbind(degred$wcon,degred$ageweight))
      
      Lix=which(contactsperson$comply==0)
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.6)
      #CTF1=0.6
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_allschool_80CT=rbind(rind_allschool_80CT,cbind(degred$wcon,degred$ageweight))
      
      CTF1=min(MAXCON/(sum(contactsperson$Numbers[Lix])/length(unique(contactsperson$PID[Lix]))),0.2)
      #CTF1=0.2
      degred=degreeperson(contactsperson,infchild,CTF1,covidsec)
      rind_allschool_20CT=rbind(rind_allschool_20CT,cbind(degred$wcon,degred$ageweight))
    }
    
    #bootstrap individual reproduction numbers and calculate confidence interval
    myboot=boot(rind,statistic=myweightedmean,R=NR)
    Rt[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_5611,statistic=myweightedmean,R=NR)
    Rt5611_noCT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt5611_noCT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_primary,statistic=myweightedmean,R=NR)
    Rtprimary_noCT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtprimary_noCT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_allschool,statistic=myweightedmean,R=NR)
    Rtall_noCT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtall_noCT[x,2:3]=mult1*bci$percent[1,4:5]
    
    ###20% CT
    myboot=boot(rind_closed_20CT,statistic=myweightedmean,R=NR)
    Rtclosed_20CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtclosed_20CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_5611_20CT,statistic=myweightedmean,R=NR)
    Rt5611_20CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt5611_20CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_primary_20CT,statistic=myweightedmean,R=NR)
    Rtprimary_20CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtprimary_20CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_allschool_20CT,statistic=myweightedmean,R=NR)
    Rtall_20CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtall_20CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    ###80% CT
    myboot=boot(rind_closed_80CT,statistic=myweightedmean,R=NR)
    Rtclosed_80CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtclosed_80CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_5611_80CT,statistic=myweightedmean,R=NR)
    Rt5611_80CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt5611_80CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_primary_80CT,statistic=myweightedmean,R=NR)
    Rtprimary_80CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtprimary_80CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    myboot=boot(rind_allschool_80CT,statistic=myweightedmean,R=NR)
    Rtall_80CT[x,1]=mult1*mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rtall_80CT[x,2:3]=mult1*bci$percent[1,4:5]
    
    x=x+1
  }
  
  
  
  # PLOTTING ---------------------------------
  
  x=100*COMP
  y2=rev(Rt[,3]); y1=rev(Rt[,2]) 
  plot(x,Rt[,1],ylim=range(Rt,0,Rt5611_noCT,Rtprimary_noCT,Rtall_noCT,3),pch=19,xlab="",col="white",
       ylab="",cex.lab=1.5,cex.axis=1.3,cex=1.0)
  grid()
  if(covidsec==0)abline(v=20,lty=2)
  abline(h=1,col="grey",lwd=3)
  
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)
  
  y2=rev(Rt5611_noCT[,3]); y1=rev(Rt5611_noCT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
  
  y2=rev(Rtprimary_noCT[,3]); y1=rev(Rtprimary_noCT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[5],border=NA)
  
  
  y2=rev(Rtall_noCT[,3]); y1=rev(Rtall_noCT[,2]) #schools partially open
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
  
  if(covidsec==1){
    mtext("0%",side=3,line=0.5,cex=1.2,bg=rgb(0.8,0.8,0.8))
    text(20,0.7,"2")
    text(97,2.8,"1")
    legend('topleft','a',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  if(covidsec==0.75){
    legend('topleft','d',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  if(covidsec==0.5){
    legend('topleft','g',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
 
  y2=rev(Rtclosed_20CT[,3]); y1=rev(Rtclosed_20CT[,2]) 
  plot(x,Rt[,1],ylim=range(Rt,0,Rt5611_noCT,Rtprimary_noCT,Rtall_noCT,3),pch=19,xlab="",col="white",
       ylab="",cex.lab=1.5,cex.axis=1.3,cex=1.3)
  

  grid()
  abline(h=1,col="grey",lwd=3)
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)
  
  y2=rev(Rt5611_20CT[,3]); y1=rev(Rt5611_20CT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
  
  y2=rev(Rtprimary_20CT[,3]); y1=rev(Rtprimary_20CT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[5],border=NA)
  
  y2=rev(Rtall_20CT[,3]); y1=rev(Rtall_20CT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
  
  if(covidsec==1){mtext("20%",side=3,line=0.5,cex=1.2,bg="grey")}
  
  if(covidsec==1){
    legend('topleft','b',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  if(covidsec==0.75){
    legend('topleft','e',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  if(covidsec==0.5){
    legend('topleft','h',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  
  x=100*COMP
  y2=rev(Rtclosed_80CT[,3]); y1=rev(Rtclosed_80CT[,2]) 
  plot(x,Rt[,1],ylim=range(Rt,0,Rt5611_noCT,Rtprimary_noCT,Rtall_noCT,3),pch=19,xlab="",col="white",
       ylab="",cex.lab=1.5,cex.axis=1.3,cex=1.3)
  #mtext("% active work & leisure contacts",side=1,line=2.5,cex=1.3)
  grid()
  abline(h=1,col="grey",lwd=3)
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)
  
  y2=rev(Rt5611_80CT[,3]); y1=rev(Rt5611_80CT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
  
  y2=rev(Rtprimary_80CT[,3]); y1=rev(Rtprimary_80CT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[5],border=NA)
  
  y2=rev(Rtall_80CT[,3]); y1=rev(Rtall_80CT[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)

  
  mtext(paste(100*(1-covidsec),"%",sep=""),side=4,line=0.5,cex=1.2,col="black",las=2)
  
  if(covidsec==1){mtext("60%",side=3,line=0.5,cex=1.2)}
  
  if(covidsec==1){
    legend('topleft','c',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  if(covidsec==0.75){
    legend('topleft','f',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
  if(covidsec==0.5){
    legend('topleft','i',cex=2,bty="n",inset=c(-0.15,-0.1))
  }
}
mtext(expression(paste("Reproduction number, ",R[t])),side=2,cex=1.5,outer=TRUE,line=1.5)
mtext("% Contacts traced",side=3,outer=TRUE,cex=1.3,line=0.5)
mtext("Covid-security",side=4,outer=TRUE,cex=1.3,line=2.2)
mtext("% active work & leisure contacts",side=1,cex=1.5,outer=TRUE,line=1.5)


####CREATE LEGEND AT BOTTOM

par(xpd=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
legend(x=-240,y=0.7,c('Schools closed',
                   '50% 5-11 year olds at school', 
                   '5-11 year olds at school',
                   'All schools open'),
       col=mycols2[c(6,2,5,4,7)],pt.bg=mycols2[c(6,2,5,4,7)],pch=22,pt.cex=2,cex=1.5,bg="white",ncol=2,y.intersp=1.2)

par(mar=rep(1.7,4),mfrow=c(3,3),oma=c(9,4,2,4),xpd=FALSE)



###END



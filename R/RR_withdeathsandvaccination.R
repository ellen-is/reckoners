#this is messy code but it is our first attempt to include vaccination into the ready reckoners
#it also has all the multiplying constants sorted out, so you can plot, for instance the distribution of individual reproduction numbers. 
#to do: comment and neaten code
#Ellen Brooks-Pollock January 2021

library(tidyverse)
library(readxl)
library(boot)

Nuk=67e6
agestats=data.frame(agegroup=c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49',
                               '50-54','55-59','60-64','65-69','70-74','75-79','80+'),
                    N=c(0.058,0.063,0.058,0.055,0.063,0.068,0.068,0.066,0.061,0.068,0.070,0.064,0.054,0.050,0.049,0.033,0.049),
                    IFR=c(0.0,0.0,0.0,0.0001,0.0002,0.0004,0.0006,0.001,0.001,0.002,0.005,0.008,0.017,0.027,0.043,0.062,0.096),
                    Hosp=c(0.001,0.001,0.001,0.002,0.005,0.01,0.016,0.023,0.029,0.039,0.058,0.072,0.102,0.117,0.146,0.177,0.18),
                    HFR=c(0.038,0.038,0.038,0.038,0.038,0.038,0.038,0.04,0.045,0.056,0.078,0.113,0.169,0.232,0.291,0.348,0.535)
)
IFR =c(0.0,0.0,0.0,0.0001,0.00015,0.0002,0.0004,0.0006,0.001,0.001,0.002,0.005,0.008,0.017,0.027,0.043,0.062,0.096,0.096,0.096)
Eng5_2021=c(3359,3716,3695,2010,1340,3612,3947,4074,3951,3739,3724,4085,4048,3501,2997,3016,2255,1545,966,434+117+15)*1000

myweightedsum=function(DD,ix)
{
  return(sum(DD[ix,1],na.rm=T))
}

myweightedmean=function(DD,ix)
{
  return(sum(DD[ix,1]^2*DD[ix,2])/sum(DD[ix,2]))
}
myweightedmean2=function(DD,ix)
{
  return(sum(DD[ix,1]*DD[ix,2])/sum(DD[ix,2]))
}

myweightedmean3=function(DD,ix)
{
  L1=length(ix)
  print("don't use weightedmean3")
  return(sum(DD[ix,1]/L1))
}

mymean=function(DD,ix)
{
  return(mean(DD[ix,1]))
}

mortality = c(rep(agestats$IFR,each=5),rep(agestats$IFR[17],15))
mortality2 = 2*mortality
hosprate = c(rep(agestats$Hosp,each=5),rep(agestats$Hosp[17],15))
hosprate2 = 2*hosprate
popsize = Nuk*c(rep(agestats$N,each=5)/5,rep((1-sum(rep(agestats$N,each=5)/5)),12))
popsize5 = Nuk*agestats$N


getsinf = function(rzero)
{
  sfinal = 0*rzero
  for(i in 1:length(rzero))
  {
    r=rzero[i]
    if(r<1){sfinal[i]=1.0}
    sinf = seq(0.001,1,0.001)
    r0=log(sinf)/(sinf-1)
    closestR=which.min(abs(r0-r))
    sfinal[i] = sinf[closestR]
  }
  return(sfinal)
}

degreeperson_calcR = function(degred,input,mult)
{  
  seroprevalence = rep(input$seroprev,max(contactsperson$age))
  
  degred %>%
    mutate(age5 = ifelse(age<=17,1+floor(age/5),ifelse(age<=90,2+floor(age/5),2+floor(90/5)))) %>%
    mutate(rd2 = recdose2[age5]) %>%
    mutate(rd1 = recdose1[age5] + input$seroprev*(1 - recdose1[age5] - rd2)) %>%
    mutate(rd0 = (1 - recdose1[age5] - rd2)*(1-input$seroprev)) %>%
    mutate(susceptibility = rd0 + rd1*(1-input$ve_inf) + rd2*(1-input$ve_inf2)) %>%
    mutate(infectiousness = rd0 + rd1*(1-input$ve_trans) + rd2*(1-input$ve_trans2)) %>%
    #mutate(infectiousness = (rd0 + rd1*(1-input$ve_inf)*(1-input$ve_trans) + rd2*(1-input$ve_trans2)*(1-input$ve_inf2))/(rd0+rd1*(1-input$ve_inf)+rd2*(1-input$ve_inf2))) %>%
    mutate(deathrate = mortality[age]*(rd0 + rd1*(1-input$ve_deathred1) + rd2*(1-input$ve_deathred2))) %>%
    mutate(deathrate2 = mortality2[age]*(rd0 + rd1*(1-input$ve_deathred1) + rd2*(1-input$ve_deathred2))) %>%
    mutate(ppp = agenumbers$ppp[age5]) %>%
    mutate(rsquare1 = mult[1]*wcon*wcon,rsquare2=mult[2]*wcon*wcon,rsquare3=mult[3]*wcon*wcon) %>%
    mutate(rind1 = mult[1]*wcon*sum(wcon),rind2=mult[2]*wcon*sum(wcon),rind3=mult[3]*wcon*sum(wcon)) %>%
    mutate(cj = wcon/sum(wcon)) %>%
    mutate(receivingrisk1 = sum(cj*(1-sigmaj1)),receivingrisk2 = sum(cj*(1-sigmaj2)),receivingrisk3 = sum(cj*(1-sigmaj3))) %>%
    mutate(sigmaj1 = exp(-rind1*receivingrisk1),sigmaj2 = exp(-rind2*receivingrisk2),sigmaj3 = exp(-rind3*receivingrisk3)) %>%
    mutate(rvsquare1 = mult[1]*susceptibility*infectiousness*wcon*wcon,
         rvsquare2=mult[2]*susceptibility*infectiousness*wcon*wcon,
         rvsquare3=mult[3]*susceptibility*infectiousness*wcon*wcon) %>%
    mutate(rpop1=sum(rvsquare1), rpop2=sum(rvsquare2), rpop3 = sum(rvsquare3)) %>%
    #mutate(rvac1 = mult[1]*wcon*susceptibility/sum(wcon),rvac2=mult[2]*wcon*susceptibility/sum(wcon),rvac3=mult[3]*wcon*susceptibility/sum(wcon)) %>%
    mutate(rvac1_unvac = mult[1]*wcon,rvac1_dose1=mult[1]*wcon*(1-input$ve_inf),rvac1_dose2=mult[1]*wcon*(1-input$ve_inf2)) %>%
    mutate(rvac2_unvac = mult[2]*wcon,rvac2_dose1=mult[2]*wcon*(1-input$ve_inf),rvac2_dose2=mult[2]*wcon*(1-input$ve_inf2)) %>%
    mutate(rvac3_unvac = mult[3]*wcon,rvac3_dose1=mult[3]*wcon*(1-input$ve_inf),rvac3_dose2=mult[3]*wcon*(1-input$ve_inf2)) %>%
    mutate(rvac_unvac = mult[3]*wcon,rvac_dose1=mult[3]*wcon*(1-input$ve_inf),rvac_dose2=mult[3]*wcon*(1-input$ve_inf2)) %>%
    mutate(cj_unvac = wcon, cj_d1 = (1-input$ve_trans)*wcon, cj_d2 = (1-input$ve_trans2)*wcon) %>%
    mutate(RR_unvac = sum(cj_unvac*(1-sigmaj_unvac)),RRvac_d1 = sum(cj_d1*(1-sigmaj_d1)),RRvac_d2 = sum(cj_d2*(1-sigmaj_d2)))%>%
    mutate(RR = sum(rd0*cj_unvac*(1-s3_unvac)+ rd1*cj_d1*(1-s3_d1) + rd2*cj_d2*(1-s3_d2)))%>%
    mutate(RR3 = sum(rd0*cj_unvac*(1-s3_unvac)+ rd1*cj_d1*(1-s3_d1) + rd2*cj_d2*(1-s3_d2)))%>%
    mutate(RR2 = sum(rd0*cj_unvac*(1-s2_unvac)+ rd1*cj_d1*(1-s2_d1) + rd2*cj_d2*(1-s2_d2)))%>%
    mutate(RR1 = sum(rd0*cj_unvac*(1-s1_unvac)+ rd1*cj_d1*(1-s1_d1) + rd2*cj_d2*(1-s1_d2)))%>%
    mutate(s3_unvac = exp(-rvac3_unvac*RR3),
           s3_d1 = exp(-rvac3_dose1*RR3),
           s3_d2 = exp(-rvac3_dose2*RR3)) %>%
    mutate(s2_unvac = exp(-rvac2_unvac*RR2),
           s2_d1 = exp(-rvac2_dose1*RR2),
           s2_d2 = exp(-rvac2_dose2*RR2)) %>%
    mutate(s1_unvac = exp(-rvac1_unvac*RR1),
           s1_d1 = exp(-rvac1_dose1*RR1),
           s1_d2 = exp(-rvac1_dose2*RR1)) %>%
    mutate(sigmaj_unvac = exp(-rvac3_unvac*RR),
           sigmaj_d1 = exp(-rvac3_dose1*RR),
           sigmaj_d2 = exp(-rvac3_dose2*RR)) %>%
    #mutate(cases3 = ppp*(rd0*(1-sigmaj_unvac) + (1-input$ve_inf)*rd1*sigmaj_d1 + (1-input$ve_inf2)*rd2*sigmaj_d2)) %>%
    #mutate(cases3 = ppp*(rd0*(1-s3_unvac) + (1-input$ve_inf)*rd1*(1-s3_d1) + (1-input$ve_inf2)*rd2*(1-s3_d2))) %>%
    #mutate(cases2 = ppp*(rd0*(1-s2_unvac) + (1-input$ve_inf)*rd1*(1-s2_d1) + (1-input$ve_inf2)*rd2*(1-s2_d2))) %>%
    #mutate(cases1 = ppp*(rd0*(1-s1_unvac) + (1-input$ve_inf)*rd1*(1-s1_d1) + (1-input$ve_inf2)*rd2*(1-s1_d2))) %>%
    #mutate(D3 = mortality2[age5]*ppp*(rd0*(1-sigmaj_unvac) + (1-input$ve_inf)*rd1*sigmaj_d1*(1-input$ve_deathred1) + (1-input$ve_inf2)*rd2*sigmaj_d2)*(1-input$ve_deathred2)) %>%
    #mutate(D3 = mortality2[age5]*ppp*(rd0*(1-s3_unvac) + (1-input$ve_inf)*rd1*(1-s3_d1)*(1-input$ve_deathred1) + (1-input$ve_inf2)*rd2*(1-s3_d2))*(1-input$ve_deathred2)) %>%
    #mutate(D2 = mortality2[age5]*ppp*(rd0*(1-s2_unvac) + (1-input$ve_inf)*rd1*(1-s2_d1)*(1-input$ve_deathred1) + (1-input$ve_inf2)*rd2*(1-s2_d2))*(1-input$ve_deathred2)) %>%
    #mutate(D1 = mortality2[age5]*ppp*(rd0*(1-s1_unvac) + (1-input$ve_inf)*rd1*(1-s1_d1)*(1-input$ve_deathred1) + (1-input$ve_inf2)*rd2*(1-s1_d2))*(1-input$ve_deathred2)) %>%
    mutate(FStest3 = (rd0*(1-s3_unvac) + rd1*(1-s3_d1) + rd2*(1-s3_d2))) %>%
    mutate(cases3 = ppp*(rd0*(1-s3_unvac) + rd1*(1-s3_d1) + rd2*(1-s3_d2))) %>%
    mutate(cases2 = ppp*(rd0*(1-s2_unvac) + rd1*(1-s2_d1) + rd2*(1-s2_d2))) %>%
    mutate(cases1 = ppp*(rd0*(1-s1_unvac) + rd1*(1-s1_d1) + rd2*(1-s1_d2))) %>%
    #mutate(D3 = mortality2[age5]*ppp*(rd0*(1-sigmaj_unvac) + (1-input$ve_inf)*rd1*sigmaj_d1*(1-input$ve_deathred1) + (1-input$ve_inf2)*rd2*sigmaj_d2)*(1-input$ve_deathred2)) %>%
    mutate(D3 = mortality2[age5]*ppp*(rd0*(1-s3_unvac) + rd1*(1-s3_d1)*(1-input$ve_deathred1) + rd2*(1-s3_d2))*(1-input$ve_deathred2)) %>%
    mutate(D2 = mortality2[age5]*ppp*(rd0*(1-s2_unvac) + rd1*(1-s2_d1)*(1-input$ve_deathred1) + rd2*(1-s2_d2))*(1-input$ve_deathred2)) %>%
    mutate(D1 = mortality2[age5]*ppp*(rd0*(1-s1_unvac) + rd1*(1-s1_d1)*(1-input$ve_deathred1) + rd2*(1-s1_d2))*(1-input$ve_deathred2)) %>%
    mutate(rvac1 = mult[1]*wcon*susceptibility,rvac2=mult[2]*wcon*susceptibility,rvac3=mult[3]*wcon*susceptibility) %>%
    mutate(cjvac = infectiousness*wcon) %>%
    mutate(RRvac1 = sum(cjvac*(1-sigmajvac1)),RRvac2 = sum(cjvac*(1-sigmajvac2)),RRvac3 = sum(cjvac*(1-sigmajvac3)))%>%
    #mutate(RRvac1 = sum(cjvac*(1-sigmajvac1)*susceptibility),RRvac2 = sum(cjvac*(1-sigmajvac2)*susceptibility),RRvac3 = sum(cjvac*(1-sigmajvac3)*susceptibility))%>%
    #mutate(RRvac1 = sum(mult[1]*cjvac*(1-sigmajvac1)*susceptibility),RRvac2 = sum(mult[2]*cjvac*(1-sigmajvac2)*susceptibility),RRvac3 = sum(mult[3]*cjvac*(1-sigmajvac3)*susceptibility))%>%
    mutate(sigmajvac1 = exp(-rvac1*RRvac1),
           sigmajvac2 = exp(-rvac2*RRvac2),
           sigmajvac3 = exp(-rvac3*RRvac3)) %>%
    mutate(sigmajvac1b = getsinf(rpop1),
           sigmajvac2b = getsinf(rpop2),
           sigmajvac3b = getsinf(rpop3)) %>%
  mutate(FS1 = (1-sigmajvac1),FS2=(1-sigmajvac2), FS3=(1-sigmajvac3)) %>%
    mutate(suspopcases = (rd0 + rd1*(1-input$ve_inf) + rd2*(1-input$ve_inf2))) %>%
    mutate(cases1b = FS1*ppp*suspopcases,cases2b=FS2*ppp*suspopcases, cases3b=FS3*ppp*suspopcases) %>%
     # mutate(cases1 = FS1*ppp*suspopcases,cases2=FS2*ppp*suspopcases) %>%
    mutate(suspop = (rd0 + rd1*(1-input$ve_inf)*(1-input$ve_deathred1) + rd2*(1-input$ve_inf2)*(1-input$ve_deathred2))) %>%
    mutate(D1b = FS1*ppp*mortality[age]*suspop,
           D2b = FS2*ppp*mortality2[age]*suspop)%>%
    mutate(D1b = FS1*ppp*mortality[age]*suspop,
           D2b = FS2*ppp*mortality2[age]*suspop, 
           D3b = FS3*ppp*mortality2[age]*suspop)%>%
    mutate(D1unvac = FS1*ppp*mortality[age]*rd0,
           D2unvac = FS2*ppp*mortality2[age]*rd0, 
           D3unvac = FS3*ppp*mortality2[age]*rd0) %>%
    mutate(D1r1 = FS1*ppp*mortality[age]*rd1*(1-input$ve_inf)*(1-input$ve_deathred1),
           D2r1 = FS2*ppp*mortality2[age]*rd1*(1-input$ve_inf)*(1-input$ve_deathred1), 
           D3r1 = FS3*ppp*mortality2[age]*rd1*(1-input$ve_inf)*(1-input$ve_deathred1)) %>%
    mutate(D1r2 = FS1*ppp*mortality[age]*rd2*(1-input$ve_inf2)*(1-input$ve_deathred2),
           D2r2 = FS2*ppp*mortality2[age]*rd2*(1-input$ve_inf2)*(1-input$ve_deathred2), 
           D3r2 = FS3*ppp*mortality2[age]*rd2*(1-input$ve_inf2)*(1-input$ve_deathred2))-> degred
  
  return(degred)
}

degreeperson = function(contactsperson,input,mult)
{
  seroprevalence = rep(input$seroprev,max(contactsperson$age))
  
  contactsperson %>%
    filter(dcat>=0)%>%
    mutate(trans=ifelse(age<=11 | Home==1,1,1-input$covidsec)) %>%
    mutate(infect = ifelse(age<=11,input$infchild,ifelse(age<=18,max(1,2*input$infchild)*trans,trans))) %>%
    mutate(Num1 = ifelse(contacttrace==TRUE,(Numbers)*(1-comply)*infect,Numbers*(1-comply)*infect)) %>%
    mutate(homenum = Num1*Home, worknum=Num1*WorkSchool, 
           othernum=Num1*OtherLeisure+Num1*Travel)%>%
    group_by(PID)%>%
    summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
              totdcat=sum(dcat),
              wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),
              ageweight=mean(ageweight))%>%
    mutate(rsquare = wcon*wcon/sum(wcon),aweight=ageweight/sum(ageweight)) %>%
    mutate(rind1 = mult[1]*wcon*sum(wcon),rind2=mult[2]*wcon*sum(wcon),rind3=mult[3]*wcon*sum(wcon))-> degred

  degred$sigmaj1=0.5; degred$sigmaj2=0.5; degred$sigmaj3=0.5; 
  degred$sigmajvac1=0.5; degred$sigmajvac2=0.5; degred$sigmajvac3=0.5; 
  degred$sigmajvac1b=0.5; degred$sigmajvac2b=0.5; degred$sigmajvac3b=0.5; 
  degred$sigmaj_unvac=0.5; degred$sigmaj_d1=0.5; degred$sigmaj_d2=0.5; 
  degred$s3_unvac=0.5; degred$s3_d1=0.5; degred$s3_d2=0.5; 
  degred$s2_unvac=0.5; degred$s2_d1=0.5; degred$s2_d2=0.5; 
  degred$s1_unvac=0.5; degred$s1_d1=0.5; degred$s1_d2=0.5; 
  degred = degreeperson_calcR(degred,input,mult)
  
  return(degred)
  
}

degreeperson0 = function(contactsperson,input)
{
  contactsperson %>%
    filter(dcat>=0)%>%
    mutate(infect = ifelse(age<=11,input$infchild,ifelse(age<=18,max(1,2*input$infchild),1))) %>%
    mutate(Num1 = ifelse(contacttrace==TRUE,(Numbers)*(1-comply)*infect,Numbers*(1-comply)*infect)) %>%
    mutate(homenum = Num1*Home, worknum=Num1*WorkSchool, 
           othernum=Num1*OtherLeisure+Num1*Travel)%>%
    group_by(PID)%>%
    summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
              totdcat=sum(dcat),
              wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),
              ageweight=mean(ageweight))%>%
    mutate(rsquare = wcon*wcon,aweight=ageweight/sum(ageweight))%>%
    mutate(rind =wcon*sum(wcon))-> degred
  
  return(degred)
}

mycols2 = rainbow(3,start=3/6,end=5/6,alpha=0.5)
mycols3 = rainbow(3,start=0/6,end=2/6,alpha=0.5)
load(file="data/uptakedose1.RData")
load(file="data/uptakedose2.RData")
Eng5_2021=c(3359,3716,3695,2010,1340,3612,3947,4074,3951,3739,3724,4085,4048,3501,2997,3016,2255,1545,966,434+117+15)*1000 #ONS 2021 population projections
IFR =c(0.0,0.0,0.0,0.0001,0.00015,0.0002,0.0004,0.0006,0.001,0.001,0.002,0.005,0.008,0.017,0.027,0.043,0.062,0.096,0.096,0.096)
mydate=as.Date("2021-05-20")
ix1=which(uptakedose1$vdate==mydate)
ix2=which(uptakedose2$vdate==mydate)
recdose1or2 = uptakedose1[ix1,3:22]
recdose2 = uptakedose2[ix2,3:22]
recdose1 = recdose1or2 - recdose2



load(file="data/contactsperson_withkids6.Rdata")
contactsperson$age[contactsperson$age==-1] = sample(contactsperson$age[contactsperson$age>0],size=sum(contactsperson$age==-1))
contactsperson %>%
  mutate(hassymp = ifelse(age<25,0.25,0.75),comply=0,contacttrace=FALSE) -> contactsperson

contactsperson %>%
  group_by(PID) %>%
  summarise(age=mean(age))%>%
  ungroup() %>%
  mutate(age5 = ifelse(age<=17,1+floor(age/5),ifelse(age<=90,2+floor(age/5),2+floor(90/5)))) %>%
  group_by(age5) %>%
  summarise(sumage=n()) %>%
  mutate(nuk = Eng5_2021[age5]) %>%
  mutate(ppp = nuk/sumage) ->agenumbers

COMP=seq(0,1,0.1)
myinterval=0.95
NR=100
#end data pre-processing

plotRt = function(Rt,Rt2,Rt3,maxY)
{
  maxY = max(maxY,Rt,Rt2,Rt3)
  x=100*COMP
  y2=rev(Rt[,3]); y1=rev(Rt[,2]) 
  
  plot(x,Rt[,1],ylim=range(0,maxY),
       pch=19,col="white",xaxt="n",yaxt="n")
  
  grid()
  abline(h=1,col="grey",lwd=3)
  #lines(x,rev(Rt[,1]))
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[1],border=NA)
  
  y2=rev(Rt2[,3]); y1=rev(Rt2[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
  #lines(x,rev(Rt2[,1]))
  
  y2=rev(Rt3[,3]); y1=rev(Rt3[,2])
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[3],border=NA)
  #lines(x,rev(Rt3[,1]))
  return(maxY)
}

plotextinction = function(Rt,Rt2,Rt3,maxY)
{
  x=100*COMP
  y1a = rev(1/Rt[,3]); y1b = rev(1/Rt[,2])
  y2a = rev(1/Rt2[,3]); y2b = rev(1/Rt2[,2])
  y3a = rev(1/Rt3[,3]); y3b = rev(1/Rt3[,2])
  
  y1a[y1a<0]=0; y1b[y1b<0]=0;
  y2a[y2a<0]=0; y2b[y2b<0]=0;
  y3a[y3a<0]=0; y3b[y3b<0]=0;
  
  y1a[y1a>1]=1; y1b[y1b>1]=1;
  y2a[y2a>1]=1; y2b[y2b>1]=1;
  y3a[y3a>1]=1; y3b[y3b>1]=1;
  
  plot(x,y3b,ylim=range(0,1),
       pch=19,col="white",xaxt="n",yaxt="n")
  
  grid()
  polygon(c(x,rev(x)),c(y3b,rev(y3a)),col=mycols2[3],border=NA)
  polygon(c(x,rev(x)),c(y2b,rev(y2a)),col=mycols2[2],border=NA)
  polygon(c(x,rev(x)),c(y1b,rev(y1a)),col=mycols2[1],border=NA)

  return(1)
}

plotfinalsize = function(Rt,Rt2,Rt3,maxY)
{
  x=100*COMP
  
  y1=rev(Rt3[,3]+1); y2=rev(Rt3[,2]+1) 
  maxY = signif(max(maxY,Rt,Rt2,Rt3),2)
  
  plot(x,y2,
       pch=19,col="white",ylim=range(1+Rt,1+Rt2,1+Rt3,maxY),xaxt="n",log="y",yaxt="n")
  
  grid()
  
  y1=rev(Rt3[,3]+1); y2=rev(Rt3[,2]+1)
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[3],border=NA)
  
  y1=rev(Rt2[,3]+1); y2=rev(Rt2[,2]+1)
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)

  y1=rev(Rt[,3]+1); y2=rev(Rt[,2]+1) 
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[1],border=NA)  
  return(maxY)
}

plotfinalsize_bystatus = function(D3U,D31,D32,maxY)
{
  x=100*COMP
  
  y0=0*D3U[,1]+1
  #y1=rev(D3U[,1]); y2=rev(D3U[,1]+D31[,1]); y3=rev(D3U[,1]+D31[,1]+D32[,1]) 
  y1=rev(D3U[,1]); y2=rev(D31[,1]); y3=rev(D32[,1]) 
  maxY = max(maxY,y0,y1,y2,y3)
  
  plot(x,y3,
       pch=19,col="white",ylim=range(1,maxY),xaxt="n",log="y",yaxt="n")
  
  grid()
  
  polygon(c(x,rev(x)),c(y0,rev(y1)),col=mycols3[3],border=NA)
  polygon(c(x,rev(x)),c(y1,rev(y2)),col=mycols3[2],border=NA)
  polygon(c(x,rev(x)),c(y2,rev(y3)),col=mycols3[1],border=NA)  

  return(maxY)
}




reckoner_variant_finalsize = function(input)
{
  baseR=c(3,5,7)
  degall=degreeperson0(contactsperson,input)
  
  mult=rep(0,3)
  mult[1]=baseR[1]/sum(degall$rsquare)
  mult[2]=baseR[2]/sum(degall$rsquare)
  mult[3]=baseR[3]/sum(degall$rsquare)
  
  #check works with boot
  #myboot=boot(degall[,c("rsquare","age")],statistic=myweightedsum,R=NR) 
  #mult[2]*mean(myboot$t)
  #mult[3]*mean(myboot$t)
  #mult[1]*mean(myboot$t)
  #perfect! 
  rind=mult[1]*degall[,c("rind")]

  rind %>%
    mutate(rind = ifelse(rind>10,10,rind)) %>%
    ggplot(aes(x=rind)) + 
    geom_histogram(colour="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept=1),
               color="black", linetype="dashed", size=1) + 
    ylab("Count") + xlab("Individual reproduction number") +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  
  ggsave("outputs/inidividualRvalues.png",width=7,height=5,units="in")
  mean(rind$rind)
  
  Rt=matrix(0,nrow=length(COMP),ncol=3)
  Rt2=matrix(0,nrow=length(COMP),ncol=3)
  Rt3=matrix(0,nrow=length(COMP),ncol=3)
  
  D1=matrix(0,nrow=length(COMP),ncol=3)
  D2=matrix(0,nrow=length(COMP),ncol=3)
  D3=matrix(0,nrow=length(COMP),ncol=3)
  
  cases1=matrix(0,nrow=length(COMP),ncol=3)
  cases2=matrix(0,nrow=length(COMP),ncol=3)
  cases3=matrix(0,nrow=length(COMP),ncol=3)
  
  D1U=matrix(0,nrow=length(COMP),ncol=3)
  D2U=matrix(0,nrow=length(COMP),ncol=3)
  D3U=matrix(0,nrow=length(COMP),ncol=3)
  
  D11=matrix(0,nrow=length(COMP),ncol=3)
  D21=matrix(0,nrow=length(COMP),ncol=3)
  D31=matrix(0,nrow=length(COMP),ncol=3)
  
  D12=matrix(0,nrow=length(COMP),ncol=3)
  D22=matrix(0,nrow=length(COMP),ncol=3)
  D32=matrix(0,nrow=length(COMP),ncol=3)
  
  #find the contact over 18
  ixO18=which(contactsperson$age>18)
  #find the school contacts
  ixU18=which(contactsperson$age<=18 & contactsperson$WorkSchool==1)
  #find leisure contacts for children
  ixU18_O=which(contactsperson$age<=18 & contactsperson$WorkSchool==0)
  #find home contacts
  ixH = which(contactsperson$Home==1)
  numreps=1
  cc=0; 
  x=1
  for(cc in COMP)
  {
    rind1=c()
    rind2=c()
    rind3=c()
    for(r in 1:numreps)
    {
      #contact tracing flag set to 0
      contactsperson$contacttrace = 0*contactsperson$dcat
      myrands = runif(length(contactsperson$hassymp),min=0,max=1)
      #contact tracing flag based on probability of symptoms
      contactsperson$contacttrace = myrands < contactsperson$hassymp

      #all schools back
      #compliance flag: 0 means contact takes place; 1 means it doesn't
      contactsperson$comply = 0*contactsperson$dcat
      contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
      contactsperson$comply[ixU18]=0 #school contacts
      contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
      contactsperson$comply[ixH]=0 #home contacts
      
      degred=degreeperson(contactsperson,input,mult)
      
      degred = degreeperson_calcR(degred,input,mult)
      for(r2 in 1:10){degred = degreeperson_calcR(degred,input,mult)}
      
      #mean(1-degred$s3_d1)
      #1-input$ve_inf
      
      #mean(1-degred$s3_d2)
      #1-input$ve_inf2
      #mean(degred$sigmajvac3)
      #mean(degred$sigmajvac3b)
      
      sum(degred$cases3)/1e6
      sum(degred$D3)
      
      rind1=rbind(rind1,degred[,c("rvsquare1","ageweight")])
      rind2=rbind(rind2,degred[,c("rvsquare2","ageweight")])
      rind3=rbind(rind3,degred[,c("rvsquare3","ageweight")])
      
      
      
    }
    myboot=boot(rind1,statistic=myweightedsum,R=NR) 
    Rt[x,1]=mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt[x,2:3]=bci$percent[1,4:5]
    
    myboot=boot(rind2,statistic=myweightedsum,R=NR) 
    Rt2[x,1]=mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt2[x,2:3]=bci$percent[1,4:5]
    
    
    myboot=boot(rind3,statistic=myweightedsum,R=NR) 
    Rt3[x,1]=mean(myboot$t)
    bci=boot.ci(myboot,type="perc",conf=myinterval)
    Rt3[x,2:3]=bci$percent[1,4:5]
  
    myboot=boot(cbind(degred[,"cases1"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    cases1[x,1]=mean(myboot$t)
    cases1[x,2:3]= range(myboot$t)
    
    myboot=boot(cbind(degred[,"cases2"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    cases2[x,1]=mean(myboot$t)
    cases2[x,2:3]= range(myboot$t)
    
    myboot=boot(cbind(degred[,"cases3"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    cases3[x,1]=mean(myboot$t)
    cases3[x,2:3]= range(myboot$t)
  
    myboot=boot(cbind(degred[,"D1"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    D1[x,1]=mean(myboot$t)
    D1[x,2:3]= range(myboot$t)
      
    myboot=boot(cbind(degred[,"D2"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    D2[x,1]=mean(myboot$t)
    D2[x,2:3]=range(myboot$t)

    myboot=boot(cbind(degred[,"D3"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
    D3[x,1]=mean(myboot$t)
    D3[x,2:3] = range(myboot$t)
      
      ####
      myboot=boot(cbind(degred[,"D1unvac"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D1U[x,1]=mean(myboot$t)
      D1U[x,2:3]= range(myboot$t)
      
      myboot=boot(cbind(degred[,"D1r1"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D11[x,1]=mean(myboot$t)
      D11[x,2:3]=range(myboot$t)
      
      myboot=boot(cbind(degred[,"D1r2"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D12[x,1]=mean(myboot$t)
      D12[x,2:3] = range(myboot$t)
      
      myboot=boot(cbind(degred[,"D2unvac"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D2U[x,1]=mean(myboot$t)
      D2U[x,2:3]= range(myboot$t)
      
      myboot=boot(cbind(degred[,"D2r1"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D21[x,1]=mean(myboot$t)
      D21[x,2:3]=range(myboot$t)
      
      myboot=boot(cbind(degred[,"D2r2"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D22[x,1]=mean(myboot$t)
      D22[x,2:3] = range(myboot$t)
      
      myboot=boot(cbind(degred[,"D3unvac"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D3U[x,1]=mean(myboot$t)
      D3U[x,2:3]= range(myboot$t)
      
      myboot=boot(cbind(degred[,"D3r1"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D31[x,1]=mean(myboot$t)
      D31[x,2:3]=range(myboot$t)
      
      myboot=boot(cbind(degred[,"D3r2"],degred[,"ageweight"]),statistic=myweightedsum,R=NR)
      D32[x,1]=mean(myboot$t)
      D32[x,2:3] = range(myboot$t)
      
    x=x+1
  }
  
  results = data.frame(R1=Rt,R2=Rt2,R3=Rt3,cases1=cases1,cases2=cases2,cases3=cases3,D1=D1,D2=D2,D3=D3,
                       D1U=D1U,D11=D11,D12=D12,D2U=D2U,D21=D21,D22=D22,D3U=D3U,D31=D31,D32=D32)
  #plotRt(Rt,Rt2,Rt3)
  #maxY=plotfinalsize(Rt,Rt2,Rt3,maxY)
  return(results)
  
}

# # 
pdf("outputs/vaccinereckoners_finalsize.pdf",height=4,width=8)
par(mar=c(0,0,0,0),mfrow=c(5,2),oma=c(7,5,5,5))
# #reckoner(input)

#three future projections
uptakedose1 = rbind(uptakedose1,uptakedose1[197,],uptakedose1[197,],uptakedose1[197,])
uptakedose1$vdate[198] = as.Date("2021-07-20")
uptakedose1$vdate[199] = as.Date("2021-08-20")
uptakedose1$vdate[200] = as.Date("2021-09-20")

uptakedose1[198,7:12] = 0.85 #vaccinate all adults over 18 with at least one dose
uptakedose1[199,7:12] = 0.85
uptakedose1[200,3:22] = uptakedose1[199,3:22]
uptakedose1[200,6] = 0.80 #vaccinate 15, 16, 17 year olds

uptakedose2 = rbind(uptakedose2,uptakedose2[176,],uptakedose2[176,],uptakedose2[176,])
uptakedose2$vdate[177] = as.Date("2021-07-20")
uptakedose2$vdate[178] = as.Date("2021-08-20")
uptakedose2$vdate[179] = as.Date("2021-09-20")

uptakedose2[177,13:15] = 0.85 #increase uptake of 2nd dose in ages 45-59
uptakedose2[178,13:15] = 0.80
uptakedose2[178,7:12] = 0.80
uptakedose2[179,3:22] = uptakedose2[178,3:22]
uptakedose2[179,6] = 0.80 #vaccinate 15, 16, 17 year olds


mydate=as.Date("2021-05-20")
#selectdates = c(as.Date("2021-03-20"),as.Date("2021-04-20"),as.Date("2021-05-20"))
selectdates = c(as.Date("2021-05-20"),as.Date("2021-06-20"),as.Date("2021-07-20"),as.Date("2021-08-20"),as.Date("2021-09-20"))
plotnum=1
maxY=0.01
#for(aget in ages)
#resultslow=c()
resultsall=c()
d1=selectdates[1]
cs=0.25
for(d1 in selectdates)
{
  for(cs in c(0,0.25))
    {
      ix1=which(uptakedose1$vdate==d1)
      ix2=which(uptakedose2$vdate==d1)
      recdose1or2 = as.numeric(uptakedose1[ix1,3:22])*1.1
      recdose2 = as.numeric(uptakedose2[ix2,3:22])*1.1
      #recdose1or2 = as.numeric(uptakedose1[ix1,3:22])
      #recdose2 = as.numeric(uptakedose2[ix2,3:22])
      recdose1 = as.numeric(recdose1or2 - recdose2)
    
    
      input = data.frame(seroprev=0.3,infchild=0.25,covidsec=cs,CTF1=0.2,
                       ve_inf=0.34,ve_trans=0.45,ve_inf2=0.73,ve_trans2=0.45,
                       ve_deathred1=0.85,ve_deathred2=0.96,maxY=4)
    
      results=reckoner_variant_finalsize(input)
      results$scenario = plotnum
      resultsall = rbind(resultsall,results)
      #resultslow = rbind(resultslow,results)
      # if(plotnum<=2){mtext(d1,side=3)}
      # if(plotnum>=9){axis(1,cex=0.5,at=seq(20,100,40))}
      # if(plotnum<=2){mtext(paste(100*cs,"%"),side=4,line=1.0)}
      # if(plotnum%%2){mtext(d1,side=3)}
      # legend('topleft',letters[plotnum],bty="n",cex=2,inset=c(-0.1,0))
      plotnum=plotnum+1
  }
}




letters=c("a","b","c","d","e","f","g","h","i","j")
mylabs=c("2021-05-20","2021-06-20","all adults\n1+ dose","all adults\n2 doses","+15-17\nyear olds")
#for(i in c(1,10,19))
for(i in c(1))
{
  par(mar=c(0,0,0,0),mfrow=c(5,2),oma=c(5,5,5,6))
  maxY=4
  maxscenarios = max(resultsall$scenario)
  for(scenario in 1:maxscenarios)
  {
    ixs=which(resultsall$scenario==scenario)
    Rt=cbind(resultsall[ixs,i],resultsall[ixs,i+1],resultsall[ixs,i+2])
    Rt2=cbind(resultsall[ixs,(i+3)],resultsall[ixs,(i+4)],resultsall[ixs,(i+5)])
    Rt3=cbind(resultsall[ixs,(i+6)],resultsall[ixs,(i+7)],resultsall[ixs,(i+8)])
    maxY=plotRt(Rt,Rt2,Rt3,maxY)
    
    if(scenario>=9){axis(1,cex=0.5,at=seq(20,100,40))}
    if(scenario==1){mtext(paste(100*0,"%"),side=3,line=1.0)}
    if(scenario==2){mtext(paste(100*0.25,"%"),side=3,line=1.0)}
    if(scenario%%2==1){axis(2,cex=0.5,at=seq(1,maxY,maxY/4))}
    if(scenario%%2==0){mtext(mylabs[scenario/2],side=4,cex=0.8,las=1,line=0.5)}
    legend('topleft',letters[scenario],bty="n",cex=2,inset=c(-0.1,0))
  }
}
mtext("COVID-security",side=3,cex=1.3,outer=TRUE,line=3)
mtext("% active work & leisure contacts",side=1,cex=1.3,outer=TRUE,line=3)
mtext(expression(paste("Reproduction number, ",R[t])),side=2,cex=1.3,outer=TRUE,line=2.5)
legend("topright",horiz = T,col=mycols2,legend=c(3,5,7),title = expression(R[0]),pch=22,pt.bg=mycols2,cex=1.0,inset=0.01,pt.cex=2)

for(i in c(10))
{
  par(mar=c(0,0,0,0),mfrow=c(5,2),oma=c(5,5,5,6))
  maxY=4
  maxscenarios = max(resultsall$scenario)
  for(scenario in 1:maxscenarios)
  {
    ixs=which(resultsall$scenario==scenario)
    Rt=cbind(resultsall[ixs,i],resultsall[ixs,i+1],resultsall[ixs,i+2])
    Rt2=cbind(resultsall[ixs,(i+3)],resultsall[ixs,(i+4)],resultsall[ixs,(i+5)])
    Rt3=cbind(resultsall[ixs,(i+6)],resultsall[ixs,(i+7)],resultsall[ixs,(i+8)])
    maxY=plotfinalsize(Rt,Rt2,Rt3,maxY)
    
    if(scenario>=9){axis(1,cex=0.5,at=seq(20,100,40))}
    if(scenario==1){mtext(paste(100*0,"%"),side=3,line=1.0)}
    if(scenario==2){mtext(paste(100*0.25,"%"),side=3,line=1.0)}
    if(scenario%%2==1){axis(2,cex=0.5,at=c(1,maxY/10000,maxY/1000,maxY/50,maxY),labels=c(1,maxY/10000,maxY/1000,maxY/50,maxY)/1e4)}
    if(scenario%%2==0){mtext(mylabs[scenario/2],side=4,cex=0.8,las=1,line=0.5)}
    legend('topleft',letters[scenario],bty="n",cex=2,inset=c(-0.1,0))
  }
}
mtext("COVID-security",side=3,cex=1.3,outer=TRUE,line=3)
mtext("% active work & leisure contacts",side=1,cex=1.3,outer=TRUE,line=3)
mtext("Total future cases x 10,000",side=2,cex=1.3,outer=TRUE,line=2.5)
legend("topright",horiz = T,col=mycols2,legend=c(3,5,7),title = expression(R[0]),pch=22,pt.bg=mycols2,cex=1.0,inset=0.01,pt.cex=2)

for(i in c(19))
{
  par(mar=c(0,0,0,0),mfrow=c(5,2),oma=c(5,5,5,6))
  maxY=4
  maxscenarios = max(resultsall$scenario)
  for(scenario in 1:maxscenarios)
  {
    ixs=which(resultsall$scenario==scenario)
    Rt=cbind(resultsall[ixs,i],resultsall[ixs,i+1],resultsall[ixs,i+2])
    Rt2=cbind(resultsall[ixs,(i+3)],resultsall[ixs,(i+4)],resultsall[ixs,(i+5)])
    Rt3=cbind(resultsall[ixs,(i+6)],resultsall[ixs,(i+7)],resultsall[ixs,(i+8)])
    maxY=plotfinalsize(Rt,Rt2,Rt3,maxY)
    
    if(scenario>=9){axis(1,cex=0.5,at=seq(20,100,40))}
    if(scenario==1){mtext(paste(100*0,"%"),side=3,line=1.0)}
    if(scenario==2){mtext(paste(100*0.25,"%"),side=3,line=1.0)}
    if(scenario%%2==1){axis(2,cex=0.5,at=c(1,10,1000,15000))}
    if(scenario%%2==0){mtext(mylabs[scenario/2],side=4,cex=0.8,las=1,line=0.5)}
    legend('topleft',letters[scenario],bty="n",cex=2,inset=c(-0.1,0))
  }
}
mtext("COVID-security",side=3,cex=1.3,outer=TRUE,line=3)
mtext("% active work & leisure contacts",side=1,cex=1.3,outer=TRUE,line=3)
mtext("Total future deaths",side=2,cex=1.3,outer=TRUE,line=2.5)
legend("topright",horiz = T,col=mycols2,legend=c(3,5,7),title = expression(R[0]),pch=22,pt.bg=mycols2,cex=1.0,inset=0.01,pt.cex=2)

for(i in c(46))
{
  par(mar=c(0,0,0,0),mfrow=c(5,2),oma=c(5,5,5,6))
  maxY=4
  maxscenarios = max(resultsall$scenario)
  for(scenario in 1:maxscenarios)
  {
    ixs=which(resultsall$scenario==scenario)
    Rt=cbind(resultsall[ixs,i],resultsall[ixs,i+1],resultsall[ixs,i+2])
    Rt2=cbind(resultsall[ixs,(i+3)],resultsall[ixs,(i+4)],resultsall[ixs,(i+5)])
    Rt3=cbind(resultsall[ixs,(i+6)],resultsall[ixs,(i+7)],resultsall[ixs,(i+8)])
    maxY=plotfinalsize_bystatus(Rt+1,Rt2+1,Rt3+1,maxY)
    
    if(scenario>=9){axis(1,cex=0.5,at=seq(20,100,40))}
    if(scenario==1){mtext(paste(100*0,"%"),side=3,line=1.0)}
    if(scenario==2){mtext(paste(100*0.25,"%"),side=3,line=1.0)}
    if(scenario%%2==1){axis(2,cex=0.5,at=c(1,10,1000,15000))}
    if(scenario%%2==0){mtext(mylabs[scenario/2],side=4,cex=0.8,las=1,line=0.5)}
    legend('topleft',letters[scenario],bty="n",cex=2,inset=c(-0.1,0))
  }
}
mtext("COVID-security",side=3,cex=1.3,outer=TRUE,line=3)
mtext("% active work & leisure contacts",side=1,cex=1.3,outer=TRUE,line=3)
mtext("Total future deaths",side=2,cex=1.3,outer=TRUE,line=2.5)
legend("topright",horiz = F,col=mycols3[3:1],legend=c("Unvaccinated","1st dose","2nd dose"),pch=22,pt.bg=mycols3[3:1],cex=0.8,inset=0.01,pt.cex=2)



for(i in c(1))
{
  par(mar=c(0,0,0,0),mfrow=c(5,2),oma=c(5,5,5,6))
  maxY=1
  maxscenarios = max(resultsall$scenario)
  for(scenario in 1:maxscenarios)
  {
    ixs=which(resultsall$scenario==scenario)
    Rt=cbind(resultsall[ixs,i],resultsall[ixs,i+1],resultsall[ixs,i+2])
    Rt2=cbind(resultsall[ixs,(i+3)],resultsall[ixs,(i+4)],resultsall[ixs,(i+5)])
    Rt3=cbind(resultsall[ixs,(i+6)],resultsall[ixs,(i+7)],resultsall[ixs,(i+8)])
    maxY=plotextinction(Rt,Rt2,Rt3,maxY)
    
    if(scenario>=9){axis(1,cex=0.5,at=seq(20,100,40))}
    if(scenario==1){mtext(paste(100*0,"%"),side=3,line=1.0)}
    if(scenario==2){mtext(paste(100*0.25,"%"),side=3,line=1.0)}
    if(scenario%%4==1){axis(2,cex=0.5,at=c(0,0.5,1))}
    if(scenario%%4==3){axis(2,cex=0.5,at=c(0.5))}
    if(scenario%%2==0){mtext(mylabs[scenario/2],side=4,cex=0.8,las=1,line=0.5)}
    legend('bottomleft',letters[scenario],bty="n",cex=2,inset=c(-0.1,0))
  }
}
mtext("COVID-security",side=3,cex=1.3,outer=TRUE,line=3)
mtext("% active work & leisure contacts",side=1,cex=1.3,outer=TRUE,line=3)
mtext("Probability of extinction",side=2,cex=1.3,outer=TRUE,line=2.5)
legend("bottomright",horiz = T,col=mycols2,legend=c(3,5,7),title = expression(R[0]),pch=22,pt.bg=mycols2,cex=1.0,inset=0.01,pt.cex=2)


agestats=data.frame(agegroup=c('0-4','5-9','10-14','15-17',"18-19",'20-24','25-29','30-34','35-39','40-44','45-49',
                               '50-54','55-59','60-64','65-69','70-74','75-79','80-84',"85-89","90+"),N=Eng5_2021,IFR=IFR,uptake=round(as.numeric(uptakedose2[179,3:22])*100))


dev.off()
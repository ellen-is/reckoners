#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(boot)


# Data pre-processing ----
myweightedmean=function(DD,ix)
{
    return(sum(DD[ix,1]^2*DD[ix,2])/sum(DD[ix,2]))
}




degreeperson = function(contactsperson,inf_child=0.25,inf_child2=0.5,CTF1,covidsec)
{
    contactsperson %>%
        filter(dcat>=0)%>%
        #mutate(trans=ifelse(Home==1,1,covidsec)) %>%
        mutate(trans=ifelse(age<=11 | Home==1,1,covidsec)) %>%
        mutate(under18 = ifelse(age<=18,1,0)) %>%
        #mutate(infect = ifelse(age<=11,inf_child,1*trans)) %>%
        mutate(infect = ifelse(age<=11,inf_child,ifelse(age<=18,inf_child2,1*trans))) %>%
        mutate(Num1 = ifelse(contacttrace==TRUE,((1-CTF1)*Numbers)*(1-comply)*infect,Numbers*(1-comply)*infect)) %>%
        mutate(homenum = Num1*Home, worknum=Num1*WorkSchool, 
               othernum=Num1*OtherLeisure+Num1*Travel)%>%
        group_by(PID)%>%
        summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
                  totdcat=sum(dcat),
                  wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),ageweight=mean(ageweight),
                  vaccinemultiplier=mean(vaccinemultiplier),aw2=mean(aw2)) -> degred
    
    return(degred)
    
}
degreeperson0 = function(contactsperson,inf_child,inf_child2=0.5,CTF1)
{
    contactsperson %>%
        filter(dcat>=0)%>%
        mutate(under18 = ifelse(age<=18,1,0)) %>%
        mutate(infect = ifelse(age<=11,inf_child,ifelse(age<=18,inf_child2,1))) %>%
        mutate(Num1 = ifelse(contacttrace==TRUE,((1-CTF1)*Numbers)*(1-comply)*infect,Numbers*(1-comply)*infect)) %>%
        mutate(homenum = Num1*Home, worknum=Num1*WorkSchool, 
               othernum=Num1*OtherLeisure+Num1*Travel)%>%
        group_by(PID)%>%
        summarise(degree=sum(homenum)+sum(worknum)+sum(othernum),
                  totdcat=sum(dcat),
                  wcon=sum(homenum*dcat2)+sum(worknum*dcat2)+sum(othernum*dcat2),age=mean(age),ageweight=mean(ageweight)) -> degred
    
    return(degred)
    
}
mycols2=c(rgb(0,0.5,1,0.5),rgb(0,0.8,0.8,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.3,0.3,0.3,0.6))

load(file="data/contactsperson_withkids6.Rdata")

contactsperson$age[contactsperson$age==-1] = sample(contactsperson$age[contactsperson$age>0],size=sum(contactsperson$age==-1))
contactsperson %>%
    mutate(hassymp = ifelse(age<25,0.75,0.25),comply=0,contacttrace=FALSE) -> contactsperson

COMP=seq(0,1,0.1)
myinterval=0.95
NR=100


#end data pre-processing 

# Define UI for application that draws a histogram
ui <- fluidPage(
    tags$head(includeHTML(("google-analytics.html"))),
    # Application title
    titlePanel("Reckoners with vaccination"),

    # Sidebar with a slider input for number of bins 
    fluidRow(
        column(2, #changed from 4 to 2 when commeted out the ve_inf and ve_trans
               h3("Vaccination parameters:"),
            selectInput('agegroup', 'Eligible age groups:', paste(seq(65,5,-5),"+",sep="")),
            #sliderInput("ve_inf",
                        "Vaccine effectiveness against infection:",
             #           min = 0,
              #          max = 1,
               #         value = 0.6),
          #  sliderInput("ve_trans",
                        "Vaccine effectiveness against transmission:",
           #             min = 0,
            #            max = 1,
             #           value = 0.5),
            sliderInput("uptake1",
                        "Vaccine uptake in low risk groups:",
                        min = 0,
                        max = 1,
                        value = 0.6),
    style="overflow-x: scroll; overflow-y: scroll"),

        # Show a plot of the generated distribution
       # mainPanel(
    column(8,
           plotOutput("distPlot")
        )
    ),
fluidRow(
    column(4,offset=0,
           h3("Transmission parameters:"),
           sliderInput("meanR0",
                       "R0:",
                       min = 0,
                       max = 5,
                       value = 3),
           sliderInput("seroprev",
                       "Seroprevalence:",
                       min = 0,
                       max = 1,
                       value = 0.25),
           #sliderInput("infchild",
              #         "Relative infectiousness of under 11s:",
            #           min = 0,
             #          max = 1,
              #         value = 0.25),
           style="overflow-x: scroll; overflow-y: scroll"),
    column(8,offset=0,
           h3(" ", style="padding:20px;"),
           sliderInput("covidsec",
                       "COVID security:",
                       min = 0,
                       max = 1,
                       value = 0),
           sliderInput("CTF1",
                       "% contacts traced:",
                       min = 0,
                       max = 100,
                       value = 10),
           style="overflow-x: scroll; overflow-y: scroll")
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        
        degall=degreeperson0(contactsperson,0.25,min(1,2*0.25),CTF1=0.0) #replaced input$infchild with 0.25
        mult1=input$meanR0/(sum((degall$wcon^2)*degall$ageweight)/sum(degall$ageweight))
        
        uptake=seq(input$uptake1,0.98,length.out=max(contactsperson$age))
        seroprevalence = rep(input$seroprev,max(contactsperson$age))

        contactsperson %>%
            mutate(vaccinated = ifelse(age>=input$agegroup,uptake[age] + (1-uptake[age])*seroprevalence[age],seroprevalence[age])) %>%
            mutate(vaccinemultiplier = sqrt(1 - 0.5*vaccinated)) %>% #0.5 was input$ve_trans
            mutate(aw2 = (1-vaccinated*0.6)*ageweight) -> contactsperson #0.6 was input$ve_inf
        
        Rt=matrix(0,nrow=length(COMP),ncol=3)
        Rt5611_noCT=matrix(0,nrow=length(COMP),ncol=3)
        Rtprimary_noCT=matrix(0,nrow=length(COMP),ncol=3)
        Rtall_noCT=matrix(0,nrow=length(COMP),ncol=3)
        
        ixO18=which(contactsperson$age>18)
        ixU18=which(contactsperson$age<=18 & contactsperson$WorkSchool==1)
        ixU18_O=which(contactsperson$age<=18 & contactsperson$WorkSchool==0)
        ixH = which(contactsperson$Home==1)
        numreps=1
        x=1
        for(cc in COMP)
        {
            rind=c()
            rind_5611=c()
            rind_primary=c()
            rind_allschool=c()
            
            for(r in 1:numreps)
            {
                contactsperson$contacttrace = 0*contactsperson$dcat
                myrands = runif(length(contactsperson$hassymp),min=0,max=1)
                contactsperson$contacttrace = myrands < contactsperson$hassymp
                
                ####now: schools closed
                contactsperson$comply = 0*contactsperson$dcat
                contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
                contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
                contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
                contactsperson$comply[ixH]=0 #home contacts
                
                
                degred=degreeperson(contactsperson,0.25,min(1,2*0.25),CTF1=input$CTF1/100,1-input$covidsec) #replaced input$infchild with 0.25
                rind=rbind(rind,cbind(degred$vaccinemultiplier*degred$wcon,degred$aw2))  
                
                ########NO CONTACT TRACING, 20% CONTACT TRACING, 60% CONTACT TRACING
                #ages 5, 6, 11 back
                contactsperson$comply = 0*contactsperson$dcat
                contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
                contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
                contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
                contactsperson$comply[ixH]=0 #home contacts
                
                ixtrans=which((contactsperson$age==11 | contactsperson$age==5 | contactsperson$age==6 | contactsperson$age==7| contactsperson$age==8) & contactsperson$WorkSchool==1)
                contactsperson$comply[ixtrans]=0
                
                
                degred=degreeperson(contactsperson,0.25,min(1,2*0.25),CTF1=input$CTF1/100,1-input$covidsec)#replaced input$infchild with 0.25
                rind_5611=rbind(rind_5611,cbind(degred$vaccinemultiplier*degred$wcon,degred$aw2))

                
                #primary schools back
                contactsperson$comply = 0*contactsperson$dcat
                contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
                contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
                contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
                contactsperson$comply[ixH]=0 #home contacts
                
                ixtrans=which((contactsperson$age<=11 & contactsperson$age>=5) & contactsperson$WorkSchool==1)
                contactsperson$comply[ixtrans]=0
                
                degred=degreeperson(contactsperson,0.25,min(1,2*0.25),CTF1=input$CTF1/100,1-input$covidsec)#replaced input$infchild with 0.25
                rind_primary=rbind(rind_primary,cbind(degred$vaccinemultiplier*degred$wcon,degred$aw2))
                
                
                #all schools back
                contactsperson$comply = 0*contactsperson$dcat
                contactsperson$comply[sample(ixO18,size=floor(cc*length(ixO18)))]=1 #other/leisure contacts for people over 18 yo
                contactsperson$comply[sample(ixU18,size=floor(0.98*length(ixU18)))]=1 #school contacts
                contactsperson$comply[sample(ixU18_O,size=floor(cc*length(ixU18_O)))]=1 #other/leisure contacts for people under 18 yo
                contactsperson$comply[ixH]=0 #home contacts
                
                ixtrans=which((contactsperson$age<=18 & contactsperson$age>=5) & contactsperson$WorkSchool==1)
                contactsperson$comply[ixtrans]=0
                
                degred=degreeperson(contactsperson,0.25,min(1,2*0.25),CTF1=input$CTF1/100,1-input$covidsec)#replaced input$infchild with 0.25
                rind_allschool=rbind(rind_allschool,cbind(degred$vaccinemultiplier*degred$wcon,degred$aw2))
                

                
            }
            
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
            
            x=x+1
        }
        x=100*COMP
        y2=rev(Rt[,3]); y1=rev(Rt[,2]) 
        par(mar=c(5,5,1,1))
        plot(x,Rt[,1],ylim=range(Rt,0,Rt5611_noCT,Rtprimary_noCT,Rtall_noCT,3),pch=19,xlab="% active work & leisure contacts",col="white",
             ylab=expression(paste("Reproduction number, ",R[t])),cex.lab=1.5,cex.axis=1.3,cex=1.0)
        grid()
        abline(h=1,col="grey",lwd=3)
        
        polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)
        
        y2=rev(Rt5611_noCT[,3]); y1=rev(Rt5611_noCT[,2])
        polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
        
        y2=rev(Rtprimary_noCT[,3]); y1=rev(Rtprimary_noCT[,2])
        polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[5],border=NA)
        
        
        y2=rev(Rtall_noCT[,3]); y1=rev(Rtall_noCT[,2]) #schools partially open
        polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
        
        #mtext(expression(paste("Reproduction number, ",R[t])),side=2,cex=1.5,outer=TRUE,line=1.5)
        #mtext("% active work & leisure contacts",side=1,cex=1.5,outer=TRUE,line=1.5)
        
        
        #par(xpd=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
        #legend('bottom',c('Schools closed',
            #                  '50% 5-11 year olds at school', 
            #                  '5-11 year olds at school',
            #                  'All schools open'),
            #   col=mycols2[c(6,2,5,4,7)],pt.bg=mycols2[c(6,2,5,4,7)],pch=22,pt.cex=2,cex=1.5,bg="white",ncol=2,y.intersp=1.2)
        
        #par(mar=rep(1.7,4),mfrow=c(3,3),oma=c(9,4,2,4),xpd=FALSE)
        
        
        
        
        #x    <- faithful[, 2]
        #bins <- seq(min(x), max(x), length.out = input$bins + 1)
        # draw the histogram with the specified number of bins
        #hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

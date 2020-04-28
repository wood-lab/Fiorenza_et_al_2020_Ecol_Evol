#### Analysis and figures for Fiorenza et al. 2020
#### "Fluid preservation causes minimal reduction of parasite detectability in fish specimens: 
#### A new approach for reconstructing parasite communities of the past?"

#Load in libraries needed
library(ggplot2)
library(lme4)
library(lmerTest)
library(metafor)
library(reshape2)
library(glmmADMB)

#Set up a standard error function
SE=function(x)
{
  sqrt(var(x)/length(x))
}

#Create a function to extract 95% confidence interval from Monte-Carlo simulations
conf95 = function(x){quantile(x,c(0.025,0.975))}


####Load in data
PollockParasites = data.frame(read.table('Data/PollockParasites.csv', header = TRUE, sep=','))
HostDat=data.frame(read.table('Data/Host Data_01.csv', header = TRUE, sep=','))
EulachonParasites = data.frame(read.table('Data/EulachonParasites.csv', header = TRUE, sep=','))
EnglishSoleParasites = data.frame(read.table('Data/EngSoleParasites.csv', header = TRUE, sep=','))
ParasiteInfo = data.frame(read.table('Data/ParasiteList.csv', header = TRUE, sep=','))

length(EulachonParasites$Host)

####Run standard linear mixed models to determine if stratified randomization achieved similar mean lengths between treatments 
polmod.length=lmer(Length~Treatment+(1|Location),data=HostDat[which(HostDat$Species=='Gadus chalcogrammus'),])
summary(polmod.length)
Eulmod.length=lmer(Length~Treatment+(1|Location),data=HostDat[which(HostDat$Species=='Thaleichthys pacificus'),])
summary(Eulmod.length)
ESmod.length=lmer(Length~Treatment+(1|Location),data=HostDat[which(HostDat$Species=='Parophrys vetulus'),])
summary(ESmod.length)


####Count up the number of fish in each species/treatment combination and give their range, mean, and SE SL
library(tidyverse)
n_full <- HostDat %>%
  group_by(Species, Treatment) %>%
  summarise(n = n(), min(Length), max(Length), mean(Length), sd(Length)/sqrt(n()))


####Fill in NAs with zeros
PollockParasites[is.na(PollockParasites)]=0
EulachonParasites[is.na(EulachonParasites)]=0
EnglishSoleParasites[is.na(EnglishSoleParasites)]=0


####Incorporate host data and parasite data together into a single data frame
PollockDat=data.frame(merge(PollockParasites,HostDat))
EulachonDat=data.frame(merge(EulachonParasites,HostDat))
EnglishSoleDat=data.frame(merge(EnglishSoleParasites,HostDat))

#how many total parasites did you find?
str(PollockDat)
numpollock<-colSums(PollockDat[,2:48])
str(EulachonDat)
numeulachon<-colSums(EulachonDat[,2:40])
str(EnglishSoleDat)
numsole<-colSums(EnglishSoleDat[,2:69])
sum(numpollock,numeulachon,numsole)

#Create lists of the parasite types found in each host species
colnames(PollockParasites)
PollockParasiteList=c('NEM.LA','NEM.LB','NEM.A','NEM.B','NEM.C','NEM.D','NEM.E','TREM.A','TREM.B','TREM.C','TREM.D','CYST.A','CYST.B','CYST.C','CYST.D','CYST.E','CYST.F','CYST.G','CYST.J','CYST.K','ACANTH.A','COPE.A','COPE.B','TRYP.A')

colnames(EulachonParasites)
EulachonParasiteList=c('NEM.LA','NEM.LB','NEM.B','NEM.C','NEM.F','CEST.B','TREM.F','TREM.G','CYST.A','CYST.E','CYST.H','CYST.I','CYST.J','ACANTH.B')

colnames(EnglishSoleParasites)
EnglishSoleParasiteList=c('CLAV','NEM.G','NEM.H','NEM.I','NEM.J','NEM.K','NEM.L1','NEM.M','NEM.LARV','NEM.ANIS','LEECH.A','CYST.A','CYST.N1','CYST.O','CYST.NEM','TREM.H','TREM.I','TREM.J','COPE.C','COPE.D','ACANTH.Z')

####Parasite data is initially separated to different organs, this code takes the sum for each parasite type across organs
for (x in 1:length(PollockParasiteList)){
  LL = PollockDat[,grep(PollockParasiteList[x], colnames(PollockDat))]
  #Create and if else choice depending on if a parasite type was found in a single organ or multiple to take the correct sum
  if (length(colnames(LL))>1){
    eval(parse(text=paste0('PollockDat$',PollockParasiteList[x],' = rowSums(LL, na.rm = TRUE)')))#Attach the row sum of the parasite type to the complete dataset
  }
  else{
    eval(parse(text=paste0('PollockDat$',PollockParasiteList[x],' = LL')))#attach the count of a parasite type  found in a single organ to the dataset
    
  }
}

for (x in 1:length(EulachonParasiteList)){
  LL = EulachonDat[,grep(EulachonParasiteList[x], colnames(EulachonDat))]
  if (length(colnames(LL))>1){
    eval(parse(text=paste0('EulachonDat$',EulachonParasiteList[x],' = rowSums(LL, na.rm = TRUE)')))
  }
  else{
    eval(parse(text=paste0('EulachonDat$',EulachonParasiteList[x],' = LL')))
    
  }
}

for (x in 1:length(EnglishSoleParasiteList)){
  LL = EnglishSoleDat[,grep(EnglishSoleParasiteList[x], colnames(EnglishSoleDat))]
  if (length(colnames(LL))>1){
    eval(parse(text=paste0('EnglishSoleDat$',EnglishSoleParasiteList[x],' = rowSums(LL, na.rm = TRUE)')))
  }
  else{
    eval(parse(text=paste0('EnglishSoleDat$',EnglishSoleParasiteList[x],' = LL')))
    
  }
}


# Determine the prevalence and abundance of parasites of each species using tidyverse

reduced_pollock<-PollockDat[,-c(2:56)]
melted_pollock<-melt(reduced_pollock,id=c("Host"))

n_prev<- melted_pollock %>%
  group_by(variable) %>%
  filter(value>0) %>%
  summarise(n = n()/109)
View(n_prev)

n_mean <- melted_pollock %>%
  group_by(variable) %>%
  summarise(n = n(), mean = mean(value), SE = sd(value)/sqrt(n()))
View(n_mean)

reduced_eulachon<-EulachonDat[,-c(2:48)]
melted_eulachon<-melt(reduced_eulachon,id=c("Host"))

n_prev<- melted_eulachon %>%
  group_by(variable) %>%
  filter(value>0) %>%
  summarise(n = n()/70)
View(n_prev)

n_mean <- melted_eulachon %>%
  group_by(variable) %>%
  summarise(n = n(), mean = mean(value), SE = sd(value)/sqrt(n()))
View(n_mean)

reduced_sole<-EnglishSoleDat[,-c(2:77)]
melted_sole<-melt(reduced_sole,id=c("Host"))

n_prev<- melted_sole %>%
  group_by(variable) %>%
  filter(value>0) %>%
  summarise(n = n()/99)
View(n_prev)

n_mean <- melted_sole %>%
  group_by(variable) %>%
  summarise(n = n(), mean = mean(value), SE = sd(value)/sqrt(n()))
View(n_mean)



#Here's another way to calculate prevalence, code written by Evan. 
PollockPrev=data.frame(Parasite=PollockParasiteList,Frozen=NA,Preserved=NA,total=NA)

for(x in 1:length(PollockParasiteList)){
  eval(parse(text=paste0('PollockPrev$Preserved[',x,']=length(PollockDat$',PollockParasiteList[x],'[which(PollockDat$',PollockParasiteList[x],'>0 & PollockDat$Treatment=="Preserved")])/54')))
  eval(parse(text=paste0('PollockPrev$Frozen[x]=length(PollockDat$',PollockParasiteList[x],'[which(PollockDat$',PollockParasiteList[x],'>0 & PollockDat$Treatment=="Fresh")])/55')))
  eval(parse(text=paste0('PollockPrev$total[x]=length(PollockDat$',PollockParasiteList[x],'[which(PollockDat$',PollockParasiteList[x],'>0)])/length(PollockDat$',PollockParasiteList[x],')')))
  
}

pol<-(PollockPrev$total>0.05)
sum(pol,na.rm=TRUE)


EulachonPrev=data.frame(Parasite=EulachonParasiteList,Frozen=NA,Preserved=NA,total=NA)

for(x in 1:length(EulachonParasiteList)){
  eval(parse(text=paste0('EulachonPrev$Preserved[',x,']=length(EulachonDat$',EulachonParasiteList[x],'[which(EulachonDat$',EulachonParasiteList[x],'>0 & EulachonDat$Treatment=="Preserved")])/27')))
  eval(parse(text=paste0('EulachonPrev$Frozen[x]=length(EulachonDat$',EulachonParasiteList[x],'[which(EulachonDat$',EulachonParasiteList[x],'>0 & EulachonDat$Treatment=="Fresh")])/43')))
  eval(parse(text=paste0('EulachonPrev$total[x]=length(EulachonDat$',EulachonParasiteList[x],'[which(EulachonDat$',EulachonParasiteList[x],'>0)])/length(EulachonDat$',EulachonParasiteList[x],')')))
  
}

eul<-(EulachonPrev$total>0.05)
sum(eul,na.rm=TRUE)


EnglishSolePrev=data.frame(Parasite=EnglishSoleParasiteList,Frozen=NA,Preserved=NA,total=NA)

for(x in 1:length(EnglishSoleParasiteList)){
  eval(parse(text=paste0('EnglishSolePrev$Preserved[',x,']=length(EnglishSoleDat$',EnglishSoleParasiteList[x],'[which(EnglishSoleDat$',EnglishSoleParasiteList[x],'>0 & EnglishSoleDat$Treatment=="Preserved")])/51')))
  eval(parse(text=paste0('EnglishSolePrev$Frozen[x]=length(EnglishSoleDat$',EnglishSoleParasiteList[x],'[which(EnglishSoleDat$',EnglishSoleParasiteList[x],'>0 & EnglishSoleDat$Treatment=="Fresh")])/48')))
  eval(parse(text=paste0('EnglishSolePrev$total[x]=length(EnglishSoleDat$',EnglishSoleParasiteList[x],'[which(EnglishSoleDat$',EnglishSoleParasiteList[x],'>0)])/length(EnglishSoleDat$',EnglishSoleParasiteList[x],')')))  
}

sole<-EnglishSolePrev$total>0.05
sum(sole,na.rm=TRUE)

#How many parasite species >0.05 in prev?
7+13+9


####For each host parastie pair, model the abundance of parasites using a negative binomial model to examine if there are significant differences between treatments
PollockModels=data.frame(ParasiteID=PollockParasiteList,Estim=NA,STD=NA,z=NA,p=NA)

for(x in c(1,4,5,8,9,12,14,21,24)){
  print(x) #Indicator of progress in the loop as well as the opportunity to know where errors arise
  eval(parse(text=paste0('',PollockParasiteList[x],'.Pollock = glmmadmb(',PollockParasiteList[x],'~Treatment+I(Length/10)+(1|Location),family="nbinom",data=PollockDat)
')))
  eval(parse(text=paste0('try=summary(',PollockParasiteList[x],'.Pollock)'))) 
  #Extract information from the models regarding the effect size, SE, Z, and p value for treatment effects
  PollockModels$Estim[x]=try$coefficients[2,1]
  PollockModels$STD[x]=try$coefficients[2,2]
  PollockModels$p[x]=try$coefficients[2,4]
  PollockModels$z[x]=try$coefficients[2,3]
  
}

EulachonModels=data.frame(ParasiteID=EulachonParasiteList,Estim=NA,STD=NA,z=NA,p=NA)

for(x in c(1,2,4,6,7,8,9)){
  print(x)  
  eval(parse(text=paste0('',EulachonParasiteList[x],'.Eulachon = glmmadmb(',EulachonParasiteList[x],'~Treatment+I(Length/10)+(1|Location),family="nbinom",data=EulachonDat)
                         ')))
  eval(parse(text=paste0('try=summary(',EulachonParasiteList[x],'.Eulachon)'))) 
  EulachonModels$Estim[x]=try$coefficients[2,1]
  EulachonModels$STD[x]=try$coefficients[2,2]
  EulachonModels$p[x]=try$coefficients[2,4]
  EulachonModels$z[x]=try$coefficients[2,3]
}

EnglishSoleModels=data.frame(ParasiteID=EnglishSoleParasiteList,Estim=NA,STD=NA,z=NA,p=NA)

for(x in c(1,2,3,4,5,7,9,11,13,14,15,16)){
  print(x)
  eval(parse(text=paste0('',EnglishSoleParasiteList[x],'.EnglishSole = glmmadmb(',EnglishSoleParasiteList[x],'~Treatment+I(Length/10)+(1|Location),family="nbinom",data=EnglishSoleDat)
                         ')))
  eval(parse(text=paste0('try=summary(',EnglishSoleParasiteList[x],'.EnglishSole)'))) 
  EnglishSoleModels$Estim[x]=try$coefficients[2,1]
  EnglishSoleModels$STD[x]=try$coefficients[2,2]
  EnglishSoleModels$p[x]=try$coefficients[2,4]
  EnglishSoleModels$z[x]=try$coefficients[2,3]
}


#Attach host identity to the extracted model information
PollockModels$Host='Pollock'
EulachonModels$Host='Eulachon'
EnglishSoleModels$Host='English Sole'
PollockModels$n=109
EulachonModels$n=70
EnglishSoleModels$n=99

#Combine the extracted model information for all hosts into a single dataframe
AllModels=rbind(PollockModels,EulachonModels,EnglishSoleModels)
#Add parasite info to the models
AllModels2=merge(AllModels,ParasiteInfo)

#how many parasite taxa did you find?
tapply(AllModels2$ParasiteID,AllModels2$Host,length)

#how many parasite taxa have a significant difference between treatments (include Bonferonni correction)?
sum(AllModels2$p<0.05/29,na.rm=TRUE)
AllModels2[38,]
AllModels2[48,]
AllModels2[50,]

#Calculate effect sizes and error for a meta-analytic model
AllMetaDat=escalc(measure='GEN', yi=Estim,sei=STD,data=AllModels2)

#Run a null meta-analytic model to see if there is a global effect of treatment
AllMetaMod.Null=rma.mv(yi=yi,V=vi,mods=~1,random = ~1|Host,data=AllMetaDat)

#Incorportate parasite information to see if group or life stage moderates the effects preservation on parsite detectability
AllMetaMod.Add=rma.mv(yi=yi,V=vi,mods=~Stage+Group,random=~1|Host,data=AllMetaDat)
AllMetaMod.Group=rma.mv(yi=yi,V=vi,mods=~Group,random=~1|Host,data=AllMetaDat)

summary(AllMetaMod.Null)
summary(AllMetaMod.Add)

####By organ meta-analysis
AllMetaDat$Pyloric=ParasiteInfo$Pyloric.Cecum[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Body=ParasiteInfo$Body[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Intestine=ParasiteInfo$Intestine[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Stomach=ParasiteInfo$Stomach[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Kidney=ParasiteInfo$Kidney[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Heart=ParasiteInfo$Heart[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Fins=ParasiteInfo$Fins[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Liver=ParasiteInfo$Liver[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Gonad=ParasiteInfo$Gonad[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Gills=ParasiteInfo$Gills[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Eye=ParasiteInfo$Eyes[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Muscle=ParasiteInfo$Muscle[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Buccal=ParasiteInfo$Buccal[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]
AllMetaDat$Skin=ParasiteInfo$Skin[match(AllMetaDat$ParasiteID,ParasiteInfo$ParasiteID)]

AllMetaMod.Organ=rma.mv(yi=yi,V=vi,mods=~Body+Pyloric+Stomach+Intestine+Kidney+Heart+Fins+Liver+Gonad+Gills+Eye+Muscle+Buccal+Skin,random=~1|Host,data=AllMetaDat)
summary(AllMetaMod.Organ)

Organ.Meta.Plot=data.frame(Organs=c('Muscle','Body','Pyloric Cecum','Stomach','Intestine','Kidney','Heart','Fins','Liver',"Gonad",'Gills','Eye','Buccal'),effect=AllMetaMod.Organ$b,LCI=AllMetaMod.Organ$ci.lb,UCI=AllMetaMod.Organ$ci.ub)



###FIGURE 1

ggplot()+
  scale_x_continuous(limits = c(0,150),name='temporal scale (years)')+
  scale_y_continuous(limits=c(-5,100),name = 'spatial scale',breaks=c(0,100),labels = c("local","global"))+
  geom_ribbon(aes(x=0:150,ymin=5,ymax=95),color=pal[1],alpha=.3,fill=pal[1],size=2)+
  geom_ribbon(aes(x=0:50,ymin=0,ymax=100),color=pal[4],alpha=.3,fill=pal[4],size=2)+
  geom_ribbon(aes(x=0:20,ymin=0,ymax=20),color=pal[7],alpha=.3,fill=pal[7],size=2)+
  geom_text(aes(x=90,y=50,label="natural history collections"),size=5)+
  geom_text(aes(x=25,y=50,label="meta-analysis"),family='sans',size=5)+
  geom_text(aes(x=10,y=11.5,label="empirical\nstudies"),size=5)+
  theme(panel.grid.major = element_blank(),axis.ticks.y=element_blank(),panel.background = element_blank(),panel.border = element_rect(fill = NA),text = element_text(family = 'sans',size=20,color="Black"))


###FIGURE 2 - See src > power_analysis_and_figure.R


###FIGURE 3
pollock.ind=PollockModels[!is.na(PollockModels$Estim),]
eulachon.ind=EulachonModels[!is.na(EulachonModels$Estim),]
sole.ind=EnglishSoleModels[!is.na(EnglishSoleModels$Estim),]  

pol.plot=ggplot(aes(x=ParasiteID,y=z),data=pollock.ind)+
  geom_rect(aes(xmin=.5,xmax=1.5,ymin=-5.5,ymax=-3.113017),fill='#ccebc5')+
  geom_rect(aes(xmin=.5,xmax=1.5,ymin=-3.113017,ymax=3.113017),fill='#ccebc5',alpha=0.1)+
  geom_rect(aes(xmin=.5,xmax=1.5,ymin=3.113017,ymax=5.5),fill='#ccebc5')+
  geom_rect(aes(xmin=1.5,xmax=2.5,ymin=-5.5,ymax=-3.113017),fill='#a8ddb5')+
  geom_rect(aes(xmin=1.5,xmax=2.5,ymin=-3.113017,ymax=3.113017),fill='#a8ddb5',alpha=0.1)+
  geom_rect(aes(xmin=1.5,xmax=2.5,ymin=3.113017,ymax=5.5),fill='#a8ddb5')+
  geom_rect(aes(xmin=2.5,xmax=5.5,ymin=-5.5,ymax=-3.113017),fill='#4eb3d3')+
  geom_rect(aes(xmin=2.5,xmax=5.5,ymin=-3.113017,ymax=3.113017),fill='#4eb3d3',alpha=0.1)+
  geom_rect(aes(xmin=2.5,xmax=5.5,ymin=3.113017,ymax=5.5),fill='#4eb3d3')+
  geom_rect(aes(xmin=5.5,xmax=9.5,ymin=-5.5,ymax=-3.113017),fill='#2b8cbe')+
  geom_rect(aes(xmin=5.5,xmax=9.5,ymin=-3.113017,ymax=3.113017),fill='#2b8cbe',alpha=0.1)+
  geom_rect(aes(xmin=5.5,xmax=9.5,ymin=3.113017,ymax=5.5),fill='#2b8cbe')+
  geom_point()+
  geom_hline(aes(yintercept=0))+
  geom_hline(aes(yintercept=-3.113017),linetype = 'dotted')+
  geom_hline(aes(yintercept=3.113017),linetype = 'dotted')+
  scale_y_continuous(limits=c(-5.5,5.5),name='',expand = c(0,0))+
  scale_x_discrete(limits=c('ACANTH.A','TRYP.A','NEM.B','NEM.C','NEM.LA','TREM.A','TREM.B','CYST.A','CYST.C'),labels=c(expression(italic("Echinorhyncus gadii")),expression(italic("Nybelinia surmenicola")),expression(paste(italic("Hysterothylacium"),' sp.')),expression(paste(italic("Contracecum"),' sp.')),expression(paste(italic("Anisakis"),' sp.')),'hemiuridean sp.','lepidapedean sp.','gill metacercariae','trematode metacercariae'),name='',expand = c(0,0))+
  labs(title=expression(italic('Gadus chalcogrammus')))+
  coord_flip()

eul.plot=ggplot(aes(x=ParasiteID,y=z),data=eulachon.ind)+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-5.5,ymax=-3.113017),fill='#a8ddb5')+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-3.113017,ymax=3.113017),fill='#a8ddb5',alpha=0.1)+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=3.113017,ymax=5.5),fill='#a8ddb5')+
  geom_rect(aes(xmin=1.5,xmax=4.5,ymin=-5.5,ymax=-3.113017),fill='#4eb3d3')+
  geom_rect(aes(xmin=1.5,xmax=4.5,ymin=-3.113017,ymax=3.113017),fill='#4eb3d3',alpha=0.1)+
  geom_rect(aes(xmin=1.5,xmax=4.5,ymin=3.113017,ymax=5.5),fill='#4eb3d3')+
  geom_rect(aes(xmin=4.5,xmax=6.5,ymin=-5.5,ymax=-3.113017),fill='#2b8cbe')+
  geom_rect(aes(xmin=4.5,xmax=6.5,ymin=-3.113017,ymax=3.113017),fill='#2b8cbe',alpha=0.1)+
  geom_rect(aes(xmin=4.5,xmax=6.5,ymin=3.113017,ymax=5.5),fill='#2b8cbe')+
  geom_point()+
  geom_hline(aes(yintercept=0))+
  geom_hline(aes(yintercept=-3.113017),linetype = 'dotted')+
  geom_hline(aes(yintercept=3.113017),linetype = 'dotted')+
  scale_y_continuous(limits=c(-5.5,5.5),name='',expand = c(0,0))+
  scale_x_discrete(limits=c('CEST.B','NEM.C','NEM.LA','NEM.LB','TREM.F','TREM.G'),labels=c('tetraphyllidean sp.',expression(paste(italic("Hysterothylacium"),' sp.')),expression(paste(italic("Anisakis"),' sp.')),expression(paste(italic("Pseudoterranova"),' sp.')),'lecithasteridean sp.',expression(paste(italic("Lecithaster"),' sp.'))),name='',expand = c(0,0))+
  labs(title=expression(italic('Thaleichthys pacificus')))+
  coord_flip()

sole.plot=ggplot(aes(x=ParasiteID,y=z),data=sole.ind)+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-5.5,ymax=-3.113017,fill='Hirundea'))+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-3.113017,ymax=3.113017,fill='Hirundea'),alpha=0.1)+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=3.113017,ymax=5.5,fill='Hirundea'))+
  geom_rect(aes(xmin=1.5,xmax=9.5,ymin=-5.5,ymax=-3.113017,fill='Nematoda'))+
  geom_rect(aes(xmin=1.5,xmax=9.5,ymin=-3.113017,ymax=3.113017,fill='Nematoda'),alpha=0.1)+
  geom_rect(aes(xmin=1.5,xmax=9.5,ymin=3.113017,ymax=5.5,fill='Nematoda'))+
  geom_rect(aes(xmin=9.5,xmax=12.5,ymin=-5.5,ymax=-3.113017,fill='Trematoda'))+
  geom_rect(aes(xmin=9.5,xmax=12.5,ymin=-3.113017,ymax=3.113017,fill='Trematoda'),alpha=0.1)+
  geom_rect(aes(xmin=9.5,xmax=12.5,ymin=3.113017,ymax=5.5,fill='Trematoda'))+
  geom_rect(aes(xmin=1,xmax=1,ymin=-10,ymax=10,fill='Acanthocephala'))+
  geom_rect(aes(xmin=1,xmax=1,ymin=-10,ymax=10,fill='Cestoda'))+
  geom_point()+
  geom_hline(aes(yintercept=0))+
  geom_hline(aes(yintercept=-3.113017),linetype = 'dotted')+
  geom_hline(aes(yintercept=3.113017),linetype = 'dotted')+
  scale_y_continuous(limits=c(-5.5,5.5),name='',expand = c(0,0))+
  scale_x_discrete(limits=c('LEECH.A','NEM.LARV','CLAV','NEM.L1','NEM.J','NEM.H','NEM.I','NEM.G','CYST.NEM','TREM.H','CYST.O','CYST.N1'),labels=c(expression(italic("Oceanobdella pallida")),'encysted larval nematode',expression(italic("Clavinema mariae")),'nematode 1','nematode 2','nematode 3','nematode 4',expression(paste(italic('Cucullanus'),' sp.')),'larval nematode',expression(paste(italic('Derogenes'),' sp.')),'fin metacercariae','trematode metacercariae'),name='',expand = c(0,0))+
  labs(title=expression(italic('Parophrys vetulus')))+
  scale_fill_manual(name='parasite group',breaks=c('Acanthocephala','Cestoda','Hirundea','Nematoda','Trematoda'),values=c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe'),labels=c('Acanthocephala','Cestoda','Hirundea','Nematoda','Trematoda'))+
  theme(legend.position = 'none')+
  coord_flip()

legend.plot=ggplot(aes(x=ParasiteID,y=z),data=sole.ind)+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-5.5,ymax=-3.113017,fill='Hirundea'))+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-3.113017,ymax=3.113017,fill='Hirundea'),alpha=0.1)+
  geom_rect(aes(xmin=0.5,xmax=1.5,ymin=3.113017,ymax=5.5,fill='Hirundea'))+
  geom_rect(aes(xmin=1.5,xmax=9.5,ymin=-5.5,ymax=-3.113017,fill='Nematoda'))+
  geom_rect(aes(xmin=1.5,xmax=9.5,ymin=-3.113017,ymax=3.113017,fill='Nematoda'),alpha=0.1)+
  geom_rect(aes(xmin=1.5,xmax=9.5,ymin=3.113017,ymax=5.5,fill='Nematoda'))+
  geom_rect(aes(xmin=9.5,xmax=12.5,ymin=-5.5,ymax=-3.113017,fill='Trematoda'))+
  geom_rect(aes(xmin=9.5,xmax=12.5,ymin=-3.113017,ymax=3.113017,fill='Trematoda'),alpha=0.1)+
  geom_rect(aes(xmin=9.5,xmax=12.5,ymin=3.113017,ymax=5.5,fill='Trematoda'))+
  geom_rect(aes(xmin=1,xmax=1,ymin=-10,ymax=10,fill='Acanthocephala'))+
  geom_rect(aes(xmin=1,xmax=1,ymin=-10,ymax=10,fill='Cestoda'))+
  geom_point()+
  geom_hline(aes(yintercept=0))+
  geom_hline(aes(yintercept=-3.113017),linetype = 'dotted')+
  geom_hline(aes(yintercept=3.113017),linetype = 'dotted')+
  scale_y_continuous(limits=c(-5.5,5.5),name='',expand = c(0,0))+
  scale_x_discrete(limits=c('LEECH.A','NEM.LARV','CLAV','NEM.L1','NEM.J','NEM.H','NEM.I','NEM.G','CYST.NEM','TREM.H','CYST.O','CYST.N1'),labels=c(expression(italic("Oceanobdella pallida")),'encysted larval nematode',expression(italic("Clavinema mariae")),'Nematode 1','Nematode 2','Nematode 3','Nematode 4',expression(paste(italic('Cucullanus'),' sp.')),'larval nematode',expression(paste(italic('Derogenes'),' sp.')),'fin metacercariae','trematode metacercariae'),name='')+
  labs(title=expression(italic('Parophrys vetulus')))+
  scale_fill_manual(name='parasite group',breaks=c('Acanthocephala','Cestoda','Hirundea','Nematoda','Trematoda'),values=c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe'),labels=c('Acanthocephala','Cestoda','Hirundea','Nematoda','Trematoda'))+
  theme(legend.direction = 'horizontal')+
  coord_flip()

arrow.plot=ggplot()+geom_path(aes(x=c(0,1),y=c(1,1)),arrow=arrow(ends='both',type='closed'))+geom_path(aes(x=c(0.01,.99),y=c(1,1)),size=2)+theme_nothing()
ind.plot.legend=get_legend(legend.plot)
ggdraw()+
  draw_plot(pol.plot,x=0.01,y=.7,height=.3,width=.99)+
  draw_plot(eul.plot,x=0.03,y=.4,height=.3,width=.97)+
  draw_plot(sole.plot,x=0.01,y=0.1,height=.3,width=.99)+
  draw_label("parasite taxon",x=0.02,y=.5,angle=90,fontfamily ='sans',size=20)+
  draw_label("standardized regression coefficient",x=0.63,y=0.085,fontfamily ='sans',size=20)+
  draw_grob(ind.plot.legend,x=0.1,y=0,height=.1,width=.8)+
  draw_plot(arrow.plot,x=.43,y=.09,height=.06,width=.4)+
  draw_label("control",x=0.395,y=0.12,fontfamily ='sans',size=19)+
  draw_label("preservation",x=0.9,y=0.12,fontfamily ='sans',size=19)


###FIGURE 4
predict.meta=data.frame(predict(AllMetaMod.Add))
AllMetaDatNArm<-AllMetaDat[complete.cases(AllMetaDat[,1:10]),]
predict.meta$Group=AllMetaDatNArm$Group
predict.meta$Stage=AllMetaDatNArm$Stage

predict.fulleffect=data.frame(predict(AllMetaMod.Null))
predict.fulleffect$Group='All'
predict.fulleffect$Stage='All'
predict.stage=data.frame(predict(AllMetaMod.Stage))
predict.stage$Stage=AllMetaDatNArm$Stage
predict.stage$Group='All'
predict.group=data.frame(predict(AllMetaMod.Group))
predict.group$Stage='All'
predict.group$Group=AllMetaDatNArm$Group

lack.data=data.frame(pred=NA,se=NA,ci.lb=NA,ci.ub=NA,cr.lb=NA,cr.ub=NA,Group=c('Leech','Acanthocephalan'),Stage='Larval')
predict.meta3=rbind(predict.meta,predict.fulleffect,predict.group,predict.stage,lack.data)
predict.meta3$Stage=factor(predict.meta3$Stage,levels = c('Larval','Adult','All'))


library(wesanderson)

wes_palette("Zissou1")
pal<-wes_palette(name = "Zissou1", 7, type = "continuous")

overallplot<-ggplot(aes(x=Group,y=pred),data=predict.fulleffect)+
  geom_crossbar(aes(ymax=ci.ub,ymin=ci.lb,fill=Group),width=.5,position=position_dodge(width = .5))+
  geom_hline(yintercept = 0)+
  ylab('')+
  xlab('')+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("All"="all taxa\nand stages"))+
  scale_fill_manual(values="darkgray")+
  theme(text=element_text(family='sans',size=20),axis.text = element_text(size=14),legend.text = element_text(size=20))

stageplot<-ggplot(aes(x=Stage,y=pred),data=predict.stage)+
  geom_crossbar(aes(ymax=ci.ub,ymin=ci.lb,fill=Stage),width=.5,position=position_dodge(width = .5))+
  geom_hline(yintercept = 0)+
  ylab('')+
  xlab('')+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Larval"="larvae","Adult"="adults"))+
  scale_fill_manual(values=pal[1:2],name='')+
  theme(text=element_text(family='sans',size=20),axis.text = element_text(size=14),legend.text = element_text(size=20))

groupplot<-ggplot(aes(x=Group,y=pred),data=predict.group)+
  geom_crossbar(aes(ymax=ci.ub,ymin=ci.lb,fill=Group),width=.5,position=position_dodge(width = .5))+
  geom_hline(yintercept = 0)+
  ylab('')+
  xlab('')+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.position = "none")+
  scale_x_discrete(limits=c("Leech","Acanthocephalan","Nematode","Cestode","Trematode"),
                   labels=c("Leech"="Hirudinea","Acanthocephalan"="Acanthocephala",
                            "Nematode"="Nematoda","Cestode"="Cestoda","Trematode"="Trematoda"))+
  scale_fill_manual(limits=c("Leech","Acanthocephalan","Nematode","Cestode","Trematode"),
                    values=c(pal[7],pal[6],pal[5],pal[4],pal[3]))+
  theme(text=element_text(family='sans',size=20),axis.text = element_text(size=14),legend.text = element_text(size=20))


library(cowplot)
ggdraw(plot=NULL,xlim=c(0,10),ylim=c(0,15))+
  draw_plot(overallplot,x=0.35,y=11.5,width=9.65,height=3)+
  draw_plot(stageplot,x=0.65,y=7.25,width=9.35,height=5)+
  draw_plot(groupplot,x=0,y=0,width=10,height=8)+
  draw_label("(a)",x=9.75,y=13.65,size=20)+
  draw_label("(b)",x=9.75,y=11.65,size=20)+
  draw_label("(c)",x=9.75,y=7.4,size=20)+
  draw_label("effect size",x=6,y=0.5,size=20)


### FIGURE 5

ggplot(aes(x=Group,y=pred),data=predict.meta3)+geom_crossbar(aes(ymax=ci.ub,ymin=ci.lb,fill=Stage),width=.5,position=position_dodge(width = .5))+geom_hline(yintercept = 0)+ylab('Effect Size')+xlab('')+coord_flip()+theme_bw()+theme(panel.grid = element_blank())+theme(legend.position = c(.8,.8))+scale_color_manual(values = c('black','black','black'))+scale_fill_manual(breaks=c('All',"Adult",'Larval'),values=c('#00e000','#67e6a1','#00e6cf'),name='Developmental Stage')+theme(text=element_text(family='sans',size=20),axis.text = element_text(size=20),legend.text = element_text(size=20))

wes_palette("Zissou1")
pal<-wes_palette(name = "Zissou1", 13, type = "continuous")

ggplot(aes(x=Organs,y=effect),data=Organ.Meta.Plot)+
  geom_crossbar(aes(ymax=UCI,ymin=LCI,fill=Organs),width=.5)+
  geom_hline(yintercept = 0)+ylab('effect size')+
  xlab('')+
  coord_flip()+
  theme_bw()+
  scale_x_discrete(limits=c(),labels=c("Muscle"="muscle","Body"="body cavity","Pyloric Cecum" = "pyloric cecae","Stomach" = "stomach",
                                       "Intestine" = "intestine","Kidney" = "kidney","Heart" ="heart","Fins" = "fins",
                                       "Liver" = "liver","Gonad" = "gonad","Gills" = "gills","Eye" = "eye",
                                       "Buccal" = "buccal cavity"))+
  scale_fill_manual(values=pal)+
  theme(legend.position = "none",text = element_text(size=20))


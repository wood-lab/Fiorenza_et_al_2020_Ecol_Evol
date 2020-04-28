#### Power analysis for Fiorenza et al. 2020
#### "Fluid preservation causes minimal reduction of parasite detectability in fish specimens: 
#### A new approach for reconstructing parasite communities of the past?"

library(lme4)
library(MASS)
library(ggplot2)
difference=seq(from=0, to=4,by=0.1)
samplesize=c(25,50,100,250,500)
total_length<-length(difference)*length(samplesize)
Poweranalysis=data.frame(Difference=rep(NA,total_length),SampleSize=rep(NA,total_length),PositiveTreatment=rep(NA,total_length))
q=seq(from=1,to=total_length,by=1)

q=1
for (a in 1:length(difference)){
  for(b in 1:length(samplesize)){
    Results=data.frame(treateffect=rep(NA,1000))
    mu1=2+difference[a]
    mu2=2
    ss=samplesize[b]
    for (x in 1:1000){
      Dat1=data.frame(treat='Frozen',count=rnbinom(n=ss,size=2,mu=mu1))
      Dat2=data.frame(treat='Formalin',count=rnbinom(n=ss,size=2,mu=mu2))
      Dat=rbind(Dat1,Dat2)
      mod.treat=glm.nb(count~treat,data=Dat)
      mod.null=glm.nb(count~1,data=Dat)
      AIC(mod.treat,mod.null)
      if (AIC(mod.treat)<AIC(mod.null)-2){
        Results$treateffect[x]=1
      }
      else{
        Results$treateffect[x]=0
      }
    }
    Poweranalysis[q,1]=difference[a]
    Poweranalysis[q,2]=samplesize[b]
    Poweranalysis[q,3]=sum(Results$treateffect)/1000
    q=q+1
    
  }
}
Poweranalysis$EffectSize=((Poweranalysis$Difference)/sqrt((((Poweranalysis$SampleSize-1)*(2+Poweranalysis$Difference+((2+Poweranalysis$Difference)^2)/2))+(Poweranalysis$SampleSize-1)*(2+((2)^2)/2))/(2*Poweranalysis$SampleSize-2)))*(1-(3/((4*Poweranalysis$SampleSize*2)-9)))

Poweranalysis$EffectSize=(Poweranalysis$Difference)/sqrt(((2+Poweranalysis$Difference+(2+Poweranalysis$Difference)^2/2)+(2+(2^2/2)))/2)


####FIGURE 2

install.packages("wesanderson")
library(wesanderson)

wes_palette("Zissou1")
pal<-wes_palette(name = "Zissou1", 5, type = "continuous")
pal<-rev(pal)

ggplot(aes(y=PositiveTreatment,x=EffectSize),data=Poweranalysis)+
  geom_point(aes(color=as.factor(SampleSize)))+
  scale_color_manual(name ="sample size\n(per treatment)", values=pal)+
  geom_hline(aes(yintercept=0.8))+
  geom_vline(aes(xintercept=0.2),linetype="dashed", col="darkgrey")+
  geom_vline(aes(xintercept=0.5),linetype="dashed", col="darkgrey")+
  geom_vline(aes(xintercept=0.8),linetype="dashed", col="darkgrey")+
  theme_minimal()+
  ylab("power")+
  xlab("effect size (Cohen's d)")+
  annotate("text", label="small",x=0.22, y=0.03, hjust=0, col="darkgrey")+
  annotate("text", label="moderate",x=0.52, y=0.03, hjust=0, col="darkgrey")+
  annotate("text", label="large",x=0.82, y=0.03, hjust=0, col="darkgrey")


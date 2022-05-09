################################################
#Analysis of phenotypic changes due to 
#hurricane disturbance
#Programmer: Miguel Acevedo
#################################################

library(ggplot2)
library(sjPlot)
library(sjmisc)
library(viridis)
library(cowplot)
library(extrafont)
library(lme4)
library(RLRsim)
library(pbkrtest)
library(effects)
library(tidyr)
library(tidyverse)
library(plyr)
library(dplyr)

#######################################################################
#Upload data
########################################################################

data <- read.csv("limbs_gundlachi_hurricane_10.csv",header=TRUE)

#########################################
#Prepare data for lizard trait analysis
#########################################
data.l <- data %>% select(uid,date,Sample,svl,sex,f1,f2,weight) %>% 
  filter(Sample == "Winter 2017" | Sample == "Winter 2018" | Sample == "Winter 2019")
data.l$Sample=revalue(data.l$Sample,c("Winter 2017" = "2017 (pre-hurricane)",
                                          "Winter 2018" = "2018 (4 months after)",
                                          "Winter 2019" = "2018 (15 months after)"))
data.l$Sample <- factor(data.l$Sample,levels=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"))

data.l$sex=revalue(data.l$sex,c("m" = "Males",
                                      "f" = "Females"))

#########################################
#Prepare data for DBH analysis
#########################################

#filter analysis for DBH analysis
#NOTE: We had no DBH data for Winter 2017 sampling, but we have data for the summer just before the hurricane.
#Therefore, we use the DBH dat afor the summer of 2017 as the pre-hurricane data for the DBH analysis.

data.dbh <- data %>% select(dbh,Sample) %>% 
  filter(Sample == "Summer 2017" | Sample == "Winter 2018" | Sample == "Winter 2019")
data.dbh$Sample=revalue(data.dbh$Sample,c("Summer 2017" = "2017 (pre-hurricane)",
                                          "Winter 2018" = "2018 (4 months after)",
                                          "Winter 2019" = "2018 (15 months after)"))
data.dbh$Sample <- factor(data.dbh$Sample,levels=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"))

######################################################################################
#SVL COMPARISON
######################################################################################
#Fixed effects model
mod.svl.fm <- lm(log(svl)~Sample+sex,data=data.l)
summary(mod.svl.fm)
tab_model(mod.svl.fm,show.se=TRUE,show.stat=TRUE)

#comparison with model without "Sample"
mod.svl.n0 <-  lm(log(svl)~sex,data=data.l)
anova(mod.svl.fm,mod.svl.n0)

#Predictions
pred.svl.m=as.data.frame(predict(mod.svl.fm,type='response',se.fit=TRUE,
                                 newdata=data.frame(Sample=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"),
                                                    sex="Males"
                                 )))
pred.svl.f=as.data.frame(predict(mod.svl.fm,type='response',se.fit=TRUE,
                                 newdata=data.frame(Sample=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"),
                                                    sex="Females"
                                 )))
pred.svl=rbind(pred.svl.m,pred.svl.f)

pred.svl$lb <- pred.svl$fit-pred.svl$se.fit*1.96
pred.svl$ub <- pred.svl$fit+pred.svl$se.fit*1.96
pred.svl$Sex=c(rep("Males",3),rep("Females",3))
pred.svl$Sample=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)",
                  "2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)")
pred.svl$Sample <- factor(pred.svl$Sample,levels=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"))

# plot
psvl <- ggplot(pred.svl, aes(Sample, exp(fit))) +
  geom_errorbar(
    aes(ymin = exp(lb), ymax = exp(ub), color = Sample),
    position = position_dodge(0.0), width = 0.2
  )+
  geom_point(aes(color = Sample,shape=Sex), position = position_dodge(0.0),size=5) +
  geom_jitter(data=data.l,aes(x=Sample,y=svl,color=Sample,shape=sex),alpha=0.30)+
  scale_color_manual(values = c( "#000000","#0072B2","#D55E00"))+ylab("snout-to-vent length (mm)")+theme_bw()+
  theme(text=element_text(family="Times New Roman"),legend.position = "top",axis.text=element_text(size=12),
        axis.title=element_text(size=12))+ 
  geom_vline(aes(xintercept=1.75), colour="black", linetype="dashed")+xlab("time")+guides(col="none")


ggsave(psvl, filename = "Fig_1.pdf", device = cairo_pdf, 
     width = 7, height = 6, units = "in")

ggsave(psvl, filename = "Fig_1.jpg", 
       width = 7, height = 6, units = "in")

ggsave(psvl, filename = "Fig_1.tiff", 
       width = 7, height = 6, units = "in",dpi=300)

##############################################################################
#Limbs analyses
##############################################################################
#f1######
mod.f1 <- lm(log(f1)~log(svl)+Sample+sex,data=data.l) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f1.re <- lmer(log(f1)~log(svl)+sex+Sample+(1|Sample:date),data=data.l, REML=FALSE)

#LRT to assess if RE are needed
exactLRT(mod.f1.re,mod.f1)

#"Significance" testing through Kenward-roger method and parametric bootstrap to confirm result
data2=data.l[complete.cases(data.l[,5:6]),] #remove rows with nas, required by PBmodcomp

mod.f1.re <- lmer(log(f1)~log(svl)+sex+Sample+(1|Sample:date),data=data2)
mod.f1.re.2 <- lmer(log(f1)~log(svl)+sex+(1|date),data=data2)

#Kenward-Roger method
kr.f1 <-KRmodcomp(mod.f1.re,mod.f1.re.2) 
summary(kr.f1)
pb.f1 <- PBmodcomp(mod.f1.re,mod.f1.re.2,nsim=5000)
pb.f1

#Make inference
summary(mod.f1.re)
tab_model(mod.f1.re,show.se=TRUE,show.stat=TRUE)

#f2###############
mod.f2 <- lm(log(f2)~log(svl)+Sample+sex,data=data.l) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f2.re <- lmer(log(f2)~log(svl)+sex+Sample+(1|Sample:date),data=data.l, REML=FALSE)

#LRT to assess if RE are needed
exactLRT(mod.f2.re,mod.f2)

mod.f2.re <- lmer(log(f2)~log(svl)+sex+Sample+(1|Sample:date),data=data2)
mod.f2.re.2 <- lmer(log(f2)~log(svl)+sex+(1|date),data=data2)

#Kenward-Roger method
kr.f2 <-KRmodcomp(mod.f2.re,mod.f2.re.2) 
summary(kr.f2)
pb.f2 <- PBmodcomp(mod.f2.re,mod.f2.re.2,nsim=5000)
pb.f2

#Make inference
summary(mod.f2.re)
tab_model(mod.f2.re,show.se=TRUE,show.stat=TRUE)


###############################
#Make limbs figure
###############################
#f1
ee1 <- Effect(c("Sample","svl","sex"),mod.f1.re)
ee1=data.frame(ee1)
ee1$Sample=factor(ee1$Sample,levels=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"))
ee1$sex=factor(ee1$sex,levels=c("Males","Females"))
ee1=ee1[-(7:15),] #remove out of range estimates for svl in females

data.l$sex=factor(data.l$sex,levels=c("Males","Females"))

plotf1 <- ggplot(ee1,
                 aes(svl,exp(fit),colour=Sample,fill=Sample))+
  facet_wrap(~sex,scales="free")+
  geom_line(aes(linetype=Sample))+
  geom_point(data=data.l,aes(x=svl,y=f1,color=Sample,shape=Sample),alpha=0.35)+
  ## colour=NA suppresses edges of the ribbon
  geom_ribbon(colour=NA,alpha=0.1,
              aes(ymin=exp(lower),ymax=exp(upper)))+
  scale_color_manual(values = c("#000000","#0072B2","#D55E00"))+
  scale_fill_manual(values = c("#000000","#0072B2","#D55E00"))+
  theme_bw()+theme(text=element_text(family="Times New Roman"),legend.position="",axis.text=element_text(size=12),
                   axis.title=element_text(size=12))+xlab("snout-to-vent length(mm)")+ylab("radius/ulna (mm)")

#f2
ee2 <- Effect(c("Sample","svl","sex"),mod.f2.re)
ee2=data.frame(ee2)
ee2$Sample=factor(ee2$Sample,levels=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"))
ee2$sex=factor(ee2$sex,levels=c("Males","Females"))
ee2=ee2[-(7:15),] #remove unrealistic estimates for svl in females

plotf2 <- ggplot(ee2,
                 aes(svl,exp(fit),colour=Sample,fill=Sample))+facet_wrap(~sex,scales="free")+
  geom_line(aes(linetype=Sample))+
  geom_point(data=data.l,aes(x=svl,y=f1,color=Sample,shape=Sample),alpha=0.35)+
  ## colour=NA suppresses edges of the ribbon
  geom_ribbon(colour=NA,alpha=0.1,
              aes(ymin=exp(lower),ymax=exp(upper)))+
  scale_color_manual(values = c("#000000","#0072B2","#D55E00"))+
  scale_fill_manual(values = c("#000000","#0072B2","#D55E00"))+
  theme_bw()+theme(text=element_text(family="Times New Roman"),legend.position="top",axis.text=element_text(size=12),
                   axis.title=element_text(size=12))+xlab("")+ylab("humerus (mm)")

pp=plot_grid(plotf2,plotf1,ncol=1, labels=c("(a)","(b)"))

ggsave(pp, filename = "Fig_2.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

ggsave(pp, filename = "Fig_2.jpg", 
       width = 6, height = 6, units = "in")

ggsave(pp, filename = "Fig_2.tiff", 
       width = 6, height = 6, units = "in",dpi=300,device=grDevices::tiff)


###########################################################################################
#DBH analysis
###########################################################################################

mod.dbh <- lm(dbh~Sample,data=data.dbh)
summary(mod.dbh)
tab_model(mod.dbh,show.se=TRUE,show.stat=TRUE)
plot_model(mod.dbh,type="pred",terms=c("Sample"))

ee1 <- Effect(c("Sample"),mod.dbh)
ee1=data.frame(ee1)
ee1$Sample <- factor(ee1$Sample,levels=c("2017 (pre-hurricane)","2018 (4 months after)","2018 (15 months after)"))

pdbh <- ggplot(ee1, aes(Sample,fit)) +
  geom_errorbar(
    aes(ymin = lower, ymax =upper),
    position = position_dodge(0.0), width = 0.2
  )+
  geom_point(aes(shape=Sample),size=5) +
  geom_jitter(data=data.dbh,aes(x=Sample,y=dbh,shape=Sample),alpha=0.30)+
  ylab("DBH (cm)")+theme_bw()+
  theme(text=element_text(family="Helvetica"),legend.position = "top",axis.text=element_text(size=12),
        axis.title=element_text(size=12))+ 
  geom_vline(aes(xintercept=1.75), colour="black", linetype="dashed")+xlab("Time")

ggsave(pdbh, filename = "plot_dbh.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

#Histogram 
pdbh.hist=ggplot(data.dbh, aes(x=dbh, fill=Sample, color=Sample)) +
  geom_histogram(position="identity",alpha=0.5)+
  scale_color_manual(values = c("#000000","#0072B2","#D55E00"))+
  scale_fill_manual(values = c("#000000","#0072B2","#D55E00"))+theme_bw()+xlab("DBH (cm)")+ylab("Frequency")
ggsave(pdbh.hist, filename = "hist_dbh.pdf", device = cairo_pdf, 
       width = 8, height = 6, units = "in")

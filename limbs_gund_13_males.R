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

data.lm <- data.l[data.l$sex=="Males",]

#Fixed effects model
mod.svl.fm <- lm(log(svl)~Sample,data=data.lm)
summary(mod.svl.fm)
tab_model(mod.svl.fm,show.se=TRUE,show.stat=TRUE)



#############################################################################
#f1######
mod.f1 <- lm(log(f1)~log(svl)+Sample,data=data.lm) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f1.re <- lmer(log(f1)~log(svl)+Sample+(1|Sample:date),data=data.lm, REML=FALSE)

#LRT to assess if RE are needed
exactLRT(mod.f1.re,mod.f1)

#"Significance" testing through Kenward-roger method and parametric bootstrap to confirm result
data2=data.lm[complete.cases(data.lm[,5:6]),] #remove rows with nas, required by PBmodcomp

mod.f1.re <- lmer(log(f1)~log(svl)+Sample+(1|Sample:date),data=data2)
mod.f1.re.2 <- lmer(log(f1)~log(svl)+(1|date),data=data2)

#Kenward-Roger method
kr.f1 <-KRmodcomp(mod.f1.re,mod.f1.re.2) 
summary(kr.f1)
pb.f1 <- PBmodcomp(mod.f1.re,mod.f1.re.2,nsim=5000)
pb.f1

#Make inference
summary(mod.f1.re)
tab_model(mod.f1.re,show.se=TRUE,show.stat=TRUE)

#f2###############
mod.f2 <- lm(log(f2)~log(svl)+Sample,data=data.lm) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f2.re <- lmer(log(f2)~log(svl)+Sample+(1|Sample:date),data=data.lm, REML=FALSE)

#LRT to assess if RE are needed
exactLRT(mod.f2.re,mod.f2)

mod.f2.re <- lmer(log(f2)~log(svl)+Sample+(1|Sample:date),data=data2)
mod.f2.re.2 <- lmer(log(f2)~log(svl)+(1|date),data=data2)

#Kenward-Roger method
kr.f2 <-KRmodcomp(mod.f2.re,mod.f2.re.2) 
summary(kr.f2)
pb.f2 <- PBmodcomp(mod.f2.re,mod.f2.re.2,nsim=5000)
pb.f2

#Make inference
summary(mod.f2.re)
tab_model(mod.f2.re,show.se=TRUE,show.stat=TRUE)

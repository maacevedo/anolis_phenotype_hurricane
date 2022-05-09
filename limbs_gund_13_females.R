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

data.lf <- data.l[data.l$sex=="Females",]
######################################################################################
#SVL COMPARISON
######################################################################################
#Fixed effects model
mod.svl.fm <- lm(log(svl)~Sample,data=data.lf)
summary(mod.svl.fm)
tab_model(mod.svl.fm,show.se=TRUE,show.stat=TRUE)


#############################################################################
#f1######
mod.f1 <- lm(log(f1)~log(svl)+Sample,data=data.lf) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f1.re <- lmer(log(f1)~log(svl)+Sample+(1|Sample:date),data=data.lf, REML=FALSE)

#Make inference
summary(mod.f1)
tab_model(mod.f1,show.se=TRUE,show.stat=TRUE)

#f2###############
mod.f2 <- lm(log(f2)~log(svl)+Sample,data=data.lf) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f2.re <- lmer(log(f2)~log(svl)+Sample+(1|Sample:date),data=data.lf, REML=FALSE)


#Make inference
summary(mod.f2)
tab_model(mod.f2,show.se=TRUE,show.stat=TRUE)

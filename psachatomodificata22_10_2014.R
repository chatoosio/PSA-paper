################################################################################### 
# R code to produce the analysis and plots of the PSA analysis of Demersal Mediterranean Stocks. Almost final version September 2014 , Author Chato Osio
#
#####################################################################################

library(stringr); library(ggplot2); library(doBy); library(gamm4);library(lme4); library(plyr); library(dplyr)

#tab<-read.csv("S:/a4a/a4a 1stmtg/gfcm_datappormethods.csv", sep=";")

setwd("S:/PSA/PSA_Final Sept 2014")

# Load PSA scores

psa <- read.csv("S:/PSA/PSA_V3/psa.csv",sep=";")


psa$DQsus <- factor(psa$DQsus, levels = c("high", "medium", "low"))
psa$DQpro <- factor(psa$DQpro, levels = c("high", "medium", "low"))
psa$DQvul <- factor(psa$DQvul, levels = c("high", "medium", "low"))

# species that are assessed
assess <-psa[psa$assess==T,]

# species with F estimation
assessz <-psa[psa$assessz==T,]

# Cartilagineous fishes
condro <-psa[psa$condro==T,]

# Species disappeared in Tyrrenian sea
disappeared <-psa[psa$disappeared==T,]

# Rename psa
#names(psa)<-c("Stock","RecordNumber","productivity","prodquality","susceptibility","susquality","vulnerability","family","supfam","vulnquality", "DQvul",  "DQpro",	"DQsus",	"IUCN",	"codici",	"vulnclass","assess", "assessz","condro","disappeared")

colnames(psa)[colnames(psa)=="pro"] <- "productivity"
colnames(psa)[colnames(psa)=="sus"] <- "susceptibility"
colnames(psa)[colnames(psa)=="specie"] <- "Stock"
colnames(psa)[colnames(psa)=="vulne"] <- "vulnerability"

############################################################################################
# decribe the data 

table(psa$vulnerability)
table(psa$productivity)
table(psa$susceptibility)

table(psa$Stock, psa$vulnerability)



####################################################################################
# TEST 3 hypothesis:
# 1 the PSA scores are consistent with the instantaneous rate of decline derived for Adriatic Elasmo

# 2 the PSA scores are consistent with the instantaneous rate of decline derived for Thyrrehenian Elasmo

# 3 Test if there is good correspondence between PSA scores and exploitation level (F/Fmsy) derived from data rich assessments

# 4 If point 3 holds, make predicitons for unassessed stocks
# What would be main assumptions and limitations of such predicitions?


# General issue: 
# 1 Reconcile  spatial scales: the PSA is not subarea specific such as Adriatic, is the comparison still valid, perhaps yes, growth parameters are similar, no evidence that these are really separate stocks

# 2 From a temporal scale, the PSA scores are from the last 4-5 years (really?), the the instantaneous yearly rate of decline over the past 50, is comparision appropriate? Susceptibility scores are from recent data (F, economic value) + selectivity is based on mesh size from recent gears. Survival should not change much over time. Vertical overlap based on trawler depth, spatial overlap atemporal. 

# 3 How do we validate the results? 

# 4 Is targeting biasing the results or since the demeral fishery is multispecies, the F levels applied to hake can be considered consistent for other stocks?

# 5 Can we incorporate data quality from PSA scores in the predicted Fratios (offset(data quality)? or model it as a variable?)

######################################################################################################
# MODEL FITTING
######################################################################################################

# Part 1 RATE OF CHANGE ~ PSA SCORES

#instMedits.csv contiene le instantaneus rates of change per Medits 1994-2005. In pratica sono i valori che corrispondono alla zona piu' ispessita delle raindrops riportati nella figura del profile likelihood (i.e i maximum likelihood estimates)

medinst<-read.csv("~/PSA/instMedits.csv")
medinst$specie<-c("Squalus blainville", "Etmopterus spinax" , "Chimaera monstrosa", "Squalus acanthias", "Mustelus asterias", "Galeus melastomus", "Scyliorhinus canicula", "Raja asterias","Raja montagui","Mustelus mustelus","Raja miraletus","Torpedo marmorata","Scyliorhinus stellaris", "Raja clavata","Myliobatis aquila","Dasyatis pastinaca","Dalatias licha","Dipturus oxyrinchus")

medinst$ratechange<-exp(medinst$betay)-1


ggplot(medinst[medinst$ratechange<10,], aes(y=ratechange,x=species))+
  geom_point()+labs(name="Year")+
  coord_flip()+
  geom_hline( yintercept=0,linetype=4 ) +
  ylab("Instantaneous decline rate")+ ggtitle ("Rate of decline in  MEDITS 1994-2005")

# Merge with PSA scores
temp0<-merge(psa, medinst, by.x="Stock", by.y="specie", all.y=TRUE) # fails to incoprorate group level estimates like Sharks Rays etc

# Fit
ggplot(temp0[temp0$ratechange<10,], aes(ratechange,susceptibility))+
  geom_point(alpha=0.5)+
  stat_smooth()+
  xlab("Instantaneous decline rate")

ggplot(temp0[temp0$ratechange<10,], aes(ratechange, productivity))+
  geom_point(alpha=0.5)+
  stat_smooth()+
  xlab("Instantaneous decline rate")

ggplot(temp0[temp0$ratechange<10,], aes(ratechange, vulnerability))+
  geom_point(alpha=0.5)+
  stat_smooth()+
  xlab("Instantaneous decline rate")



#HvarMedits.csv sono il long-term estimates per ottenere la diminuzione o aumento in percentuale devi fare exp(estimate)-1(quelli che corrispondono alla figura 5)
hvarmed<-read.csv("~/PSA/HvarMedits.csv")
# hvarmed$ratechange<-exp(hvarmed$estimate)-1 # Following F Ferretti's advice

colnames(hvarmed)[colnames(hvarmed)=="estimate"] <- "ratechange"

ggplot(hvarmed, aes(y=ratechange,x=species))+
  geom_point()+labs(name="Year")+
  coord_flip()+
  geom_hline( yintercept=0,linetype=4 ) +
  ylab("Instantaneous decline rate")+ ggtitle ("Rate of decline in HVAR-MEDITS")

# Galeus melastomus and Raja asterias appear as strong outliers in terms of delcine rate, need to decide is include or not

summary(medinst)

# Merge psa scorse with historical estimates of rate of decline in Hvar vs MEDITS

temp<-merge(psa, hvarmed, by.x="Stock", by.y="species", all.y=TRUE) # fails to incoprorate group level estimates like Sharks Rays etc
#merge(psa, hvarmed, by.x="Stock", by.y="species", all.y=TRUE)
#"Rajiformes","Squaliformes" + Elasmobranchs

ggplot(temp[temp$ratechange<10,], aes(ratechange,vulnerability))+geom_point(alpha=0.5)+stat_smooth(method="lm")
ggplot(temp[temp$ratechange<10,], aes(ratechange,susceptibility))+geom_point(alpha=0.5)+stat_smooth(method="lm")
ggplot(temp[temp$ratechange<10,], aes(ratechange,productivity))+geom_point(alpha=0.5)+stat_smooth(method="lm")



############################################################################

# Part 2 Testing PSA scores ~ Rate of Decline in Thyrreninan Sea

#coastal waters.csv 
coastal<-read.csv("~/PSA/coastal2.csv", sep=",")
#str_sub(coastal$species, start = -1, end = -1) <- ""; coastal$species
temp4<-merge(psa, coastal, by.x="Stock", by.y="species", all.y=TRUE )
temp4<- temp4[temp4$Parameter=="year",]

# temp4$ratechange<-exp(as.numeric(as.character(temp4$Estimate)))-1
#exp(medinst$betay)-1

plot(temp4$vulnerability, temp4$ratechange)
ggplot(temp4[temp4$ratechange >-0.4,], aes(ratechange,vulnerability))+geom_point()+stat_smooth(method="lm")



#deep waters.csv 
deep<-read.csv("~/PSA/deep2.csv", sep=",")
#str_sub(deep$species, start = -1, end = -1) <- ""; deep$species
temp5<-merge(psa, deep, by.x="Stock", by.y="species", all.y=TRUE)

temp5<- temp5[temp5$Parameter=="year",]

### levo la Raja polystigma
temp5<- temp5[-10,]



# Check consistency of rates of decline in Adriatic and Tyrrhenian
# If ok try to merge the last 3 dataset and retry fitting the models
# These are simple glm's, perhaps more advanced methods would fit much differently??
# Is a problem modelling a rate?  

#######################################################################################################
# merge all the rates of decline
hvarmed$zone<-factor("adriatic")
coastal$zone<-factor("tyrrhenian")
deep$zone<-factor("tyrrhenian")
# medinst$zone<-factor("adriatic")


coastal<- coastal[coastal$Parameter=="year",]
deep<- deep[deep$Parameter=="year",]
# coastal$ratechange<-exp(as.numeric(as.character(coastal$Estimate)))-1
# deep$ratechange<-exp(as.numeric(as.character(deep$Estimate)))-1



# prova<-hvarmed[,c(1,3,6,7)]
names(hvarmed)
prova<-hvarmed[,as.character(c("species","ratechange","se","zone"))]
prova2<-deep[,as.character(c("species","Estimate","StdEr","zone"))]
prova3<-coastal[,as.character(c("species","Estimate","StdEr","zone"))]
names(prova)<-c("species","Estimate", "StdEr",     "zone"   )
prova4<-rbind(prova,prova2,prova3)


bigone<-merge(psa, prova4, by.x="Stock", by.y="species", all.y=TRUE)
bigone $ Stock <- bigone $ Stock[drop=TRUE]
bigone2 <- bigone
bigone2 $ wts <- 1/bigone2 $ StdEr^2
bigone2 <- na.omit(bigone2)
bigone2 <- subset(bigone2, is.finite(wts))
bigone2$Estimate<-as.numeric(as.character(bigone2$Estimate))

# Model fitting PSA scores ~ Rate of Change

ggplot(bigone2, aes(x=susceptibility, y=Estimate))+geom_point()+stat_smooth(method="lm")+facet_grid(~zone)

# remove outliers
bigone2<-bigone2[-32,]
bigone2<-bigone2[-21,]
bigone2<-bigone2[-11,]
bigone2<-bigone2[-14,]

# fit a first round of gamm's with smoother of interaction productivity and susceptibility
mod1 <-gamm4(Estimate ~ t2(productivity, susceptibility, k=4)+Stock,
             random=~ (1|zone), weights = bigone2$wts, data=bigone2, REML=FALSE)
mod2 <-gamm4(Estimate ~ t2(productivity, susceptibility, k=4),
             random=~ (1|zone)+ (1|Stock), weights= bigone2$wts, data=bigone2, REML=FALSE)
mod3 <-gamm4(Estimate ~ t2(productivity , susceptibility, k=4)+zone,
             random=~ (1|Stock), weights=bigone2$wts,data=bigone2, REML=FALSE)
mod4 <-gamm4(Estimate ~ t2(productivity , susceptibility, k=4)+zone,
             weights=bigone2$wts,data=bigone2, REML=FALSE)

mod41 <-gamm4(Estimate ~ s(productivity, k=4)+ s(susceptibility, k=4),
              weights=bigone2$wts,data=bigone2, REML=FALSE)
mod41a <-gamm4(Estimate ~ s(productivity, k=4)+ s(susceptibility, k=4),
               random=~ (1|zone), weights=bigone2$wts,data=bigone2, REML=FALSE)

mod42 <-gamm4(Estimate ~ productivity * susceptibility,
              random=~ (1|zone), weights=bigone2$wts,data=bigone2, REML=FALSE)
mod43 <-gamm4(Estimate ~ productivity,
              random=~ (1|zone),   weights=bigone2$wts,data=bigone2, REML=FALSE)
mod44 <-gamm4(Estimate ~ susceptibility,
              random=~ (1|zone), weights=bigone2$wts,data=bigone2 , REML=FALSE)
#mod45 <-gamm4(Estimate ~ susceptibility+zone ,           weights=wts,data=bigone2)

mod31 <-gamm4(Estimate ~ s(productivity,k=4)+ s(susceptibility, k=4),
              weights=bigone2$wts,data=bigone2, REML=FALSE)

mod5 <-gam(Estimate ~ susceptibility+zone,
           weights=bigone2$wts, data=bigone2, REML=FALSE)

mod6 <-gamm4(Estimate ~ t2(productivity, susceptibility, k=4)+zone, 
             random=~ (1|Stock), weights = bigone2$wts, data=bigone2, REML=FALSE)

# Get AIC
AIC(mod1$mer,mod2$mer,mod3$mer, mod4$mer,mod31$mer, mod41$mer,mod41a$mer, mod42$mer, 
    mod43$mer, mod44$mer, mod5, mod6$mer)

# calculate Delta AIC
(AIC(mod1$mer,mod2$mer,mod3$mer, mod4$mer,mod31$mer, mod41$mer,mod41a$mer, mod42$mer, 
     mod43$mer, mod44$mer, mod5)$AIC)-   
  min(AIC(mod1$mer,mod2$mer,mod3$mer, mod4$mer,mod31$mer, mod41$mer,mod41a$mer, 
          mod42$mer, mod43$mer, mod44$mer, mod5)$AIC)
#####################################
# Alternative from Colin
models1 <- list(mod1,mod2,mod3, mod4,mod31, mod41,mod41a, mod42, 
    mod43, mod44,  mod6)
modnames1 <- c("mod1","mod2","mod3", "mod4","mod31", "mod41","mod41a", "mod42", 
    "mod43", "mod44", "mod6")

# get AIC, df, edf, and deltaAIC
tab1 <- data.frame(
  model = modnames1,
  df = sapply(models1, function(x) attr(logLik(x $ mer), "df")),
  AIC = sapply(models1, function(x) AIC(x $ mer)),
  dAIC = NA,
  edf = sapply(models, function(x) sum(x $ gam $ edf))
)
tab1 $ dAIC <- tab1 $ AIC - min(tab1 $ AIC)
tab1



#####################################
gam.check(mod5)

summary(mod5)

# Get some diagnostics
png(file="mod5_Diagnostics2.png",  width = 780, height = 780,  antialias = "cleartype")

gam.check(mod5)

dev.off()

# Extract vector of Effective Degrees of Freedom from model fits
c(summary(mod1$gam)$edf[1],
  c(summary(mod2$gam)$edf[1],summary(mod3$gam)$edf[1],summary(mod4$gam)$edf[1],
    summary(mod31$gam)$edf[1],summary(mod41$gam)$edf[1],summary(mod41a$gam)$edf[1],
    summary(mod42$gam)$edf[1],summary(mod43$gam)$edf[1],summary(mod44$gam)$edf[1])) 



##########################################################################################
# Part 3 Fitting Vulnerability scores ~ F/Fmsy

#Import all the data for current assessment results
#stock1 <- read.csv("~/PSA/STECF_summary_allAssessments2008_2012_12022013new.csv",sep=";")
stock1 <- read.csv("~/PSA/stock1.csv", sep=";")
stock1<-subset(stock1, select=-c(X, X.1, part))
#duplicated assessment for Hake in GSA11 in 2012 remove part 1,
#stock1<-stock1[-423,]
#stock1<-stock1[-507,]  #stock1[stock1$Species=="Trisopterus minutus" & stock1$GSA==9,]

stockl<-stock1
stockl$Fbar<-as.numeric(as.character(stockl$Fbar))
stockl$Fratio<-as.numeric(as.character(stockl$Fbar))/as.numeric(as.character(stockl$Fmsy))
stockl$Fratio<-as.numeric(as.character(stockl$Fratio))
stockl$year<-as.numeric(as.character(stockl$year))-1 # -1 to convert from the year of the assessment to year of biological data

#import additional groupings
species_groups_codes <- read.csv("~/PSA/species_groups_codes.csv", sep=";")


# Arrange the levels of GSAs
#stockl<-subset(stockl, !GSA %in% "29")
#gsa<- factor(c("1","5","6","7","8","9","10","11","15-16","16","17","18", "19","20","22","22-23","25") )   

#
gsa<- factor(c("1","5","6","7","8","9","10","11","12-16","15-16","16","17","18", "19","20","22","22-23","25") )   
stockl$GSA<-factor(stockl$GSA, levels=gsa)


# If new Species are addedd to the original list they need to be added here otherways will be be a level mismatch or some will be dropped 
species <- factor(c("Engraulis encrasicolus","Sardina pilchardus", "Amblyraja radiata","Aristaeomorpha foliacea",  "Aristeus antennatus","Boops boops", "Galeus melastomus","Lophius budegassa","Merlangius merlangus","Merluccius merluccius","Micromesistius poutassou", "Mullus barbatus","Mullus surmuletus", "Nephrops norvegicus","Octopus vulgaris","Pagellus erythrinus","Parapenaeus longirostris","Phycis blennoides", "Psetta maxima","Rapana venosa", "Solea solea","Spicara smaris","Spicara flexuosa", "Sprattus sprattus", "Squalus acanthias", "Squilla mantis", "Trachurus mediterraneus",  "Trisopterus minutus") ) 
stockl$Species<-factor(stockl$Species, levels=species)

stockl2<-merge(stockl, species_groups_codes, by="Species", all.x=TRUE)

# retain species for which there is at least one assessment
# stockl<-subset(stockl, Fratio>=0 & year>2008) 
stockl2<-subset(stockl2, Fratio>=0) 

# in temp2 I am bringing in pelagic and crustacean stocks, need to remove to have constant dimensions for predictions with new data set
stockl2$GROUP<-as.character(stockl2$GROUP)
stockl2<-stockl2[stockl2$GROUP=="Demersal Fish" | stockl2$GROUP=="Demersal Elasmobranch",]
summary(stockl2)
names(psa)
names(stockl)
temp2<-merge(psa, stockl2, by.x="Stock", by.y="Species", all.x=TRUE)
#names(temp2)
#temp2<-	subset(temp2, select = -c(fit, ciu, cil))

names(temp2)<-c("Stock","RecordNumber","productivity","prodquality","susceptibility","susquality","vulnerability","family","supfam","vulnquality", "DQvul",  "DQpro",  "DQsus",  "IUCN",	"codici",	"vulnclass","assess", "assessz","condro","disappeared","GSA","Common.name","year","Fbar","Fmsy","Fratio","X","code","GROUP")
head(temp2)

#plot Fratio against vulnerability
ggplot(temp2, aes(vulnerability, Fratio))+geom_point()+stat_smooth(method="lm")+facet_wrap(~GSA)

temp2$yearF<-as.factor(temp2$year)
temp2$YS<-paste(temp2$year, temp2$Stock, sep="") # creates a combination of year and stock to check if there are trends, the problem is that all the unassessed stocks have NAs and are then dropped from the model

# remove black sea stocks, area 29
#temp2$GSA <- as.character(temp2$GSA)
temp3<-temp2

#temp3$GSA<-as.factor(temp3$GSA)
summary(temp3)





## NOTE only data where there is a non-NA GSA is used
temp3 <- subset(temp3, !is.na(GSA))
temp3 $ GSA <- temp3 $ GSA[drop = TRUE]
temp3 <- temp3[temp3$GSA!="12-16",]
temp3$GSA<- factor (temp3$GSA)


# start selecting the random part

Fmod1g <-gamm4(log(Fratio) ~ t2(productivity, susceptibility, k=4)+GSA,
               random=~(1| YS) + (1|Stock),    data=temp3, REML=FALSE)

Fmod2g <-gamm4(log(Fratio) ~  t2(productivity, susceptibility, k=4)+GSA,
               random=~(1|Stock),    data=temp3, REML=FALSE)

Fmod2ag <-gamm4(log(Fratio) ~ t2(productivity, susceptibility, k=4)+GSA,
                random=~(1| YS) + (1|GSA),    data=temp3, REML=FALSE)

Fmod3g <-gamm4(log(Fratio) ~ t2(productivity,  susceptibility, k=4)+GSA,
               random=~  (1|YS),    data=temp3, REML=FALSE)

# try to drop interaction
Fmod4g <-gamm4(log(Fratio) ~ s(productivity, k=4)+ s(susceptibility,k=4)+GSA,
               random=~ (1| YS)+(1|Stock),    data=temp3, REML=FALSE, na.action=na.omit)
# drop individual terms
Fmod5g <-gamm4(log(Fratio) ~ s(susceptibility,k=4)+GSA,
               random=~ (1| YS)+(1|Stock),    data=temp3, REML=FALSE)
Fmod6g <-gamm4(log(Fratio) ~ s(productivity, k=4)+GSA, 
               random=~ (1| YS)+(1|Stock), data=temp3, REML=FALSE)
Fmod7g <-gamm4(log(Fratio) ~ s(productivity, k=4) +GSA, 
               random=~ (1|Stock), data=temp3, REML=FALSE)
Fmod8g <-gamm4(log(Fratio) ~ s(susceptibility, k=4)+GSA, 
               random=~ (1| Stock), data=temp3, REML=FALSE)
Fmod9g <-gamm4(log(Fratio) ~ s(productivity, k=4),     
               random=~ (1| YS), data=temp3, REML=FALSE)
Fmod7sg <-gamm4(log(Fratio) ~ productivity + GSA  ,
                random=~ (1|Stock), data=temp3, REML=FALSE)


models <- list(Fmod1g, Fmod2g, Fmod2ag, Fmod3g, Fmod4g,Fmod5g,Fmod6g, Fmod7g, 
               Fmod8g, Fmod9g, Fmod7sg)
modnames <- c("Fmod1g", "Fmod2g", "Fmod2ag", "Fmod3g", "Fmod4g", "Fmod5g", "Fmod6g", "Fmod7g", 
              "Fmod8g", "Fmod9g", "Fmod7sg")

# get AIC, df, edf, and deltaAIC
tab <- data.frame(
  model = modnames,
  df = sapply(models, function(x) attr(logLik(x $ mer), "df")),
  AIC = sapply(models, function(x) AIC(x $ mer)),
  dAIC = NA,
  edf = sapply(models, function(x) sum(x $ gam $ edf))
)
tab $ dAIC <- tab $ AIC - min(tab $ AIC)
tab


# refit best model with REML
Fmod7sgREML <-gamm4(log(Fratio) ~ productivity +GSA  , random=~ (1|Stock), data=temp3, REML=TRUE)

summary(Fmod7sgREML$mer)
summary(Fmod7sgREML$gam)

# Get some diagnostics
png(file="Fmod7gREML_Diagnostics.png",  width = 480, height = 480)
gam.check((Fmod7sgREML$gam))
dev.off()



# PREDICTION FUNCTION looping over GSAs to get a full dataset

predStatus <- function(newdata, model = Fmod7sgREML $ gam) {
  pred <- predict(model, newdata, se.fit = TRUE)
     
# add predictions with approximate 95% confidence interval to prediction data.frame
  newdata $ fit <- exp(pred $ fit)
  newdata $ ciu <- exp(pred $ fit + 2*pred $ se)
  newdata $ cil <- exp(pred $ fit - 2*pred $ se)

#  a <- ggplot(x, aes(productivity, fit, color=fit,shape=GSA, ymin=cil, ymax=ciu))+ geom_point(position = "jitter")+scale_colour_gradient2(midpoint=1, low="Darkgreen", high="red", name="Exploitation")+ ylab("Predicted(F/Fmsy)")+geom_errorbar(position = "jitter")+stat_smooth()+facet_grid(.~GSA)
#print(a)
  newdata
}


newdata <- psa[c("Stock", "productivity", "susceptibility")]

datafinal <- do.call(rbind, lapply(levels(temp3 $ GSA), function(gsa) cbind(newdata, GSA = gsa)))
datafinal <- datafinal[datafinal$GSA!="12-16",]# this gsa was giving problems with dimensions
datafinal$GSA<- factor (datafinal$GSA)
datafinal <- predStatus(datafinal, model = Fmod7sgREML $ gam)

# Plot predicted F/Fmsy for the unassessed stocks
ggplot(datafinal, aes(productivity, fit, color=fit, ymin=cil, ymax=ciu))+ 
    geom_point()+
    scale_colour_gradient(low="Darkgreen", high="red", name="Exploitation")+ 
    ylab("Predicted(F/Fmsy)")+xlab("Productivity")+
    geom_errorbar()+
    #stat_smooth(method = "loess")+
    facet_grid(.~GSA) +
    theme(legend.position = "bottom")+
   #scale_y_continuous(breaks=c(1,4,7,9))+
  scale_x_continuous(breaks=c(1,2,3))+
  geom_hline( yintercept=1,linetype=4 )
ggsave(last_plot(), file=paste("plotPredictFFmsy",Sys.Date(),".png"), width=14, height=8, dpi=300)

# PLot fitted against predicted for assessed stocks

datafinal2 <- predStatus(temp3, model = Fmod7sgREML $ gam)

ggplot(datafinal2, aes(productivity, fit, color = Stock, ymin = cil, ymax = ciu)) +
    ylab("Predicted(F/Fmsy)") +xlab("Productivity")+
    geom_pointrange(position = position_dodge(width = 0.9))+
    geom_point(aes(x = productivity, y = Fratio, color = Stock), shape = 4, 
               position = position_dodge(width = 0.9)) +
    facet_grid(.~GSA)+
    guides(col = guide_legend(nrow = 2)) +
    theme(legend.position = "bottom")+
  scale_x_continuous(breaks=c(1,2,3))+
  geom_hline( yintercept=1,linetype=4 )
ggsave(last_plot(), file=paste("plotPredict_vs_Observed_FFmsy",Sys.Date(),".png"), width=11, height=8, dpi=300)

# try to plot overall results by super family and area

datafinal2 <-merge(datafinal, subset(psa, select=c("Stock", "family", "supfam", "condro"), by=all))

# Plot density by GSA
ggplot(datafinal2, aes(x=fit))+geom_density(aes(fill=factor(supfam),alpha=0.5))+
  facet_grid(GSA~supfam)+ geom_vline( xintercept=1,linetype=4 ) +
  #geom_text(data=NULL, mapping=aes(x=1, y=0, label="F/Fmsy=1"), size=4, angle=90, vjust=-0.4, hjust=0)+
  xlab("Predicted F/Fmsy")+scale_x_continuous(breaks=c(1,3,5,7))+
  scale_y_continuous(breaks=c(0,1,2))+scale_fill_discrete(name="Superfamily")

# Plot density for all areas combined, allows visualization of all superfamilies
ggplot(datafinal3, aes(x=fit))+geom_density(aes(fill=(supfam),alpha=0.5))+facet_grid(GSA~., scales= "free_y")+ 
  theme(strip.text.y = element_text(size = 8, angle = 0))+
  geom_vline( xintercept=1,linetype=4 ) +
  scale_fill_discrete(name="Superfamily")+
   xlab("Predicted F/Fmsy")

ggsave(last_plot(), file=paste("plotPredictFFmsy_KernelDensity",Sys.Date(),".png"), width=9, height=8, dpi=300)




########### 

# Plot the ranking of Price*Landings+Vuln
ranking <- read.csv("S:/PSA/PSA_Final Sept 2014/ranking_land_val_vuln.csv")
names(ranking)<-c("landings", "price", "vuln", "Scientific_name", "rank")
spe <- reorder(ranking$Scientific_name,ranking$rank)
ranking <- ranking[order(ranking$rank, decreasing = TRUE),] 

# define the stocks that are assessed
temp4<-unique(subset(temp2, select= c("Stock", "assess")))
ranking2<-merge(ranking, temp4, by.x="Scientific_name", by.y="Stock")

# reorder the species based on rank
ranking2$Scientific_name <- factor(ranking2$Scientific_name, levels = ranking2$Scientific_name[order(ranking2$rank)])

ggplot(ranking2, aes(Scientific_name,rank, stat="bin"))+
  geom_bar(aes(color=assess))+coord_flip()+
scale_colour_discrete(name  = "Assessed Stock")+
  xlab("")+ylab("Rank")+
  opts(axis.text.y = theme_text(face = 'italic'))

ggsave(last_plot(), file=paste("StockRankings",Sys.Date(),".png"), width=9, height=8, dpi=300)


# Summarize results for results and discussion

head(datafinal)
xtab(~GSA+fit, datafinal)
a<-ddply(datafinal[datafinal$fit>1,], .(GSA), summarize,  n.stocks=length(unique(as.character(Stock))), .drop = TRUE)

# select stocks exploited unsustainably
a<-datafinal[datafinal$fit<1,]
a1<-merge(a, subset(psa, select=c("Stock", "family", "supfam"), by=all))
table(a1$family, a1$GSA)
table(a1$supfam[a1$GSA==20], a1$GSA[a1$GSA==20])

# explore results
aa<-datafinal[datafinal$cil<1,]
aa1<-merge(aa, subset(psa, select=c("Stock", "family", "supfam"), by=all))
table(aa1$family, aa1$GSA)

write.csv(merge(datafinal, subset(psa, select=c("Stock", "family", "supfam"), by=all)), file="predictedExploitation.csv")

datafinal[datafinal$cil>3,]
ddply(datafinal3[datafinal3$fit>1,], .(supfam,GSA), summarize,  n.stocks=length(unique(as.character(Stock))), .drop = TRUE)

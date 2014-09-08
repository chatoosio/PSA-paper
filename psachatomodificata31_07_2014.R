################################################################################### 
# R code to produce the analysis and plots of the PSA analysis of Demersal Mediterranean Stocks. Almost final version September 2014 , Author Chato Osio
#
#####################################################################################

library(stringr); library(ggplot2); library(doBy); library(gamm4);library(lme4)

#tab<-read.csv("S:/a4a/a4a 1stmtg/gfcm_datappormethods.csv", sep=";")

 setwd("S:/PSA/PSA_V3")

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
# Part 1
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
mod1 <-gamm4(Estimate ~ s(productivity * susceptibility, k=4)+Stock , random=~ (1|zone), weights = bigone2$wts, data=bigone2, REML=FALSE)
mod2 <-gamm4(Estimate ~ s(productivity * susceptibility, k=4) , random=~ (1|zone)+ (1|Stock), weights=bigone2wts,data=bigone2, REML=FALSE)
mod3 <-gamm4(Estimate ~ s(productivity * susceptibility, k=4)+zone , random=~ (1|Stock), weights=bigone2$wts,data=bigone2, REML=FALSE)
mod4 <-gamm4(Estimate ~ s(productivity * susceptibility, k=4)+zone, weights=bigone2$wts,data=bigone2, REML=FALSE)

mod41<-gamm4(Estimate ~ s(productivity, k=4)+ s(susceptibility, k=4), weights=bigone2$wts,data=bigone2, REML=FALSE)
mod41a <-gamm4(Estimate ~ s(productivity, k=4)+ s(susceptibility, k=4) , random=~ (1|zone), weights=bigone2$wts,data=bigone2, REML=FALSE)

mod42 <-gamm4(Estimate ~ productivity * susceptibility , random=~ (1|zone), weights=bigone2$wts,data=bigone2, REML=FALSE)
mod43 <-gamm4(Estimate ~ productivity  , random=~ (1|zone),   weights=bigone2$wts,data=bigone2, REML=FALSE)
mod44 <-gamm4(Estimate ~ susceptibility , random=~ (1|zone), weights=bigone2$wts,data=bigone2 , REML=FALSE)
#mod45 <-gamm4(Estimate ~ susceptibility+zone ,           weights=wts,data=bigone2)

mod31 <-gamm4(Estimate ~ s(productivity,k=4)+ s(susceptibility, k=4) ,           weights=bigone2$wts,data=bigone2, REML=FALSE)

AIC(mod1$mer,mod2$mer,mod3$mer, mod4$mer,mod31$mer, mod41$mer,mod41a$mer, mod42$mer, mod43$mer, mod44$mer)

# refit best model with REML
mod44REML <-gamm4(Estimate ~ susceptibility , random=~ (1|zone), weights=bigone2$wts,data=bigone2, REML=TRUE)
#mod44REML <-gamm4(Estimate ~ susceptibility + zone , random=~ (1|zone), weights=bigone2$wts,data=bigone2, REML=TRUE)

summary(mod44REML$mer)
summary(mod44REML$gam)

# Generate a new data frame for prediction
newdata<-data.frame(susceptibility=c(1.5, 2 , 2.5 , 3), zone=factor(c("adriatic", "tyrrhenian")))

pred <- predict(mod44REML$gam, newdata, se.fit = TRUE)

# add predictions with approximate 95% confidence interval to prediction data.frame
newdata$fit <- pred $ fit
newdata$ciu <- pred $ fit + 2*pred $ se
newdata$cil <- pred $ fit - 2*pred $ se

# we plot the data, the fit and the upper and lower confidence intervals
#plot(newdata$)
plot(bigone2$susceptibility[bigone2$zone=="adriatic"], bigone2$Estimate[bigone2$zone=="adriatic"])
lines(newdata$susceptibility, newdata$fit, col = "blue", lwd = 2)
lines(newdata$susceptibility, newdata$ciu, col = "blue", lwd = 2, lty = 2)
lines(newdata$susceptibility, newdata$cil, col = "blue", lwd = 2, lty = 2)

ggplot(newdata, aes(susceptibility, fit))+
	geom_line(aes(susceptibility, fit, fill=zone))+
	geom_ribbon(data=newdata,aes(ymin=cil,ymax=ciu, fill=zone),alpha=0.2)


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

# start selecting the random part

Fmod1g <-gamm4(Fratio ~ s(productivity, by=susceptibility, k=4)+GSA , random=~(1| YS) + (1|Stock),    data=temp3, REML=FALSE)

Fmod2g <-gamm4(Fratio ~  s(productivity, by=susceptibility, k=4)+GSA , random=~(1|Stock),    data=temp3, REML=FALSE)

#Fmod2g <-gamm4(Fratio ~ s(productivity,susceptibility, k=4)+GSA  , random=~(1| YS) + (1|GSA),    data=temp3, REML=FALSE)
# Error in checkNlevels(reTrms$flist, n = n, control) : 
#  grouping factors must have > 1 sampled level
Fmod3g <-gamm4(Fratio ~ s(productivity, by = susceptibility, k=4)+GSA  , random=~  (1|YS),    data=temp3, REML=FALSE)

anova(Fmod1g$gam, Fmod2g$gam)
anova(Fmod2g$gam, Fmod3g$gam)

# try to drop interaction
Fmod4g <-gamm4(Fratio ~ s(productivity, k=4)+ s(susceptibility,k=4)+GSA  , random=~ (1| YS)+(1|Stock),    data=temp3, REML=FALSE, na.action=na.omit)
# drop individual terms
Fmod5g <-gamm4(Fratio ~ s(susceptibility,k=4)+GSA  , random=~ (1| YS)+(1|Stock),    data=temp3, REML=FALSE)
Fmod6g <-gamm4(Fratio ~ s(productivity, k=4) +GSA  , random=~ (1| YS)+(1|Stock),    data=temp3, REML=FALSE)
Fmod7g <-gamm4(Fratio ~ s(productivity, k=4) +GSA  , random=~ (1|Stock),            data=temp3, REML=FALSE)
Fmod8g <-gamm4(Fratio ~ s(susceptibility, k=4) +GSA  , random=~ (1| Stock),              data=temp3, REML=FALSE)
Fmod9g <-gamm4(Fratio ~ s(productivity, k=4)  ,      random=~ (1| YS),              data=temp3, REML=FALSE)

summary(Fmod7g$gam)
AIC(Fmod1g$mer,Fmod3g$mer, Fmod4g$mer,Fmod5g$mer,Fmod6g$mer, Fmod7g$mer, Fmod8g$mer, Fmod9g$mer)
#aa <- as.list(formula( Fmod4g$gam), formula( Fmod5g$gam),formula( Fmod6g$gam), formula( Fmod7g$gam), formula( Fmod8g$gam), formula( Fmod9g$gam))
#do.call(formula, Fmod9g$gam)

# refit best model with REML
Fmod7gREML <-gamm4(Fratio ~ s(productivity, k=4) +GSA  , random=~ (1|Stock), data=temp3, REML=TRUE)

summary(Fmod7gREML$mer)
summary(Fmod7gREML$gam)


# PREDICTION FUNCTION looping over GSAs to get a full dataset

predStatus<-function(x, y, z){
     pred <- predict(y, x, se.fit = TRUE)
     x<-data.frame( Stock=(psa$Stock), productivity=psa$productivity, susceptibility=psa$susceptibility, GSA=factor(z)) 
     
# add predictions with approximate 95% confidence interval to prediction data.frame
    x$fit <- pred $ fit
    x$ciu <- pred $ fit + 2*pred $ se
    x$cil <- pred $ fit - 2*pred $ se

    assign(paste("area", z, sep=""), x)


      a<-ggplot(x, aes(productivity, fit, color=fit,shape=GSA, ymin=cil, ymax=ciu))+ geom_point(position = "jitter")+scale_colour_gradient2(midpoint=1, low="Darkgreen", high="red", name="Exploitation")+ ylab("Predicted(F/Fmsy)")+geom_errorbar(position = "jitter")+stat_smooth()+facet_grid(.~GSA)

  print(a)
 #   return(x)
}

predStatus(x = temp3, y = Fmod9gREML$gam, z = 9 )

# Loop over all GSA
 for ( i in stockl$GSA){
   a<-predStatus(newdata2, Fmod9gREML$gam, z= i)
   assign(paste("area", i, sep=""), a)
   return(a)
 }
datafinal<-rbind(area1, area5,area6 ,area7 ,area9 ,area10 ,area11   ,area17 ,area18 ,area19, area20, area25 )

# Combine the datasets
datafinal<-rbind(area1, area5,area6 ,area7 ,area9 ,area10 ,area11   ,area17 ,area18 ,area19, area20, area25 )



ggplot(datafinal, aes(productivity, fit, color=fit,shape=GSA, ymin=cil, ymax=ciu))+ 
    geom_point(position = "jitter")+
  scale_colour_gradient2(midpoint=3, low="Darkgreen", high="red", name="Exploitation")+ 
  ylab("Predicted(F/Fmsy)")+
  geom_errorbar(position = "jitter")+
  stat_smooth()+
  facet_grid(.~GSA)


# Plot the predicted exploitation rates by area and stock with stock names
ggplot(datafinal, aes(productivity, fit, color=fit, label=Stock))+ 
  #  geom_point(position = "jitter")+ 
  geom_text(position=position_jitter(h=1,w=1))+
  scale_colour_gradient2(midpoint=3, low="Darkgreen", high="red", name="Exploitation")+ 
  ylab("Predicted(F/Fmsy)")+
  facet_grid(.~GSA)









# Attempt prdicition
newdata2<-data.frame( Stock=(psa$Stock), productivity=psa$productivity, susceptibility=psa$susceptibility, GSA=factor(17))  

predStatus<-function(x, y){
  pred <- predict(y, x, se.fit = TRUE)

# add predictions with approximate 95% confidence interval to prediction data.frame
x$fit <- pred $ fit
x$ciu <- pred $ fit + 2*pred $ se
x$cil <- pred $ fit - 2*pred $ se

assign(paste("area", i, sep=""), x)


a<-ggplot(x, aes(productivity, fit, color=fit,shape=GSA, ymin=cil, ymax=ciu))+ geom_point(position = "jitter")+scale_colour_gradient2(midpoint=1, low="Darkgreen", high="red", name="Exploitation")+ ylab("Predicted(F/Fmsy)")+geom_errorbar(position = "jitter")+stat_smooth()+facet_grid(.~GSA)
print(a)
return(x)
}

#newdata<-rbind(newdata, newdata2)

predStatus(newdata2, Fmod8gREML$gam)



ggplot(x, aes(productivity, fit, color=fit,shape=GSA, ymin=cil, ymax=ciu))+ geom_point(position = "jitter")+scale_colour_gradient2(midpoint=3, low="Darkgreen", high="red", name="Exploitation")+ ylab("Predicted(F/Fmsy)")+
geom_errorbar(position = "jitter")+stat_smooth()+facet_grid(.~GSA)
}
#and
vis.gam(Fmod8gREML$gam, view=c("susceptibility", "productivity"), theta=35)

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	# Only working prediciton.....
	newdata2<-data.frame( Stock=(psa$Stock), productivity=psa$productivity, susceptibility=psa$susceptibility, GSA= factor(1))  

#predStatus(x = newdata2, y = Fmod9gREML$gam, z = newdata2$GSA )

for ( i in (newdata2$GSA)){
	a<-predStatus(newdata2, Fmod7gREML$gam, z = i)
	assign(paste("areaN", i, sep=""), a)
	
}

datafinal<-rbind(areaN1, areaN5,areaN6 ,areaN7 ,areaN9 ,areaN10 ,areaN11   ,areaN17 ,areaN18 ,areaN19, areaN20, areaN25)

ggplot(datafinal, aes(productivity, fit, color=fit, ymin=cil, ymax=ciu))+ 
	#  scale_colour_gradient2(midpoint=3, low="Darkgreen", high="red", name="Exploitation")+ 
	scale_colour_gradientn(colours=c("green","white","khaki2","yellow" ,"yellow2","orange","orange2","red","red2"),limits=c(0, 10),breaks=0:10, name="F/Fmsy")+
	ylab("Predicted(F/Fmsy)")+
	geom_errorbar(position = "jitter", alpha = 1/5)+
	# stat_smooth()+
	geom_hline( yintercept=1,linetype=4 )+
	geom_point(position = "jitter", size=2)+
	facet_grid(.~GSA)

ggsave(last_plot(), file=paste("plotPredictFFmsy",Sys.Date(),".png"), width=14, height=8, dpi=300)


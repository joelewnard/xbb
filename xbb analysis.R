setwd("~/Library/CloudStorage/GoogleDrive-jlewnard@berkeley.edu/My Drive/CAP KP study/xbb/data analysis")
load('completeDat.Rdata')

cur = completeDat[[1]]
quantile(cur$tVax1[which(cur$totVax=='1'&cur$sgtf=='1')],c(0.5,0.25,0.75));
quantile(cur$tVax1[which(cur$totVax=='1'&cur$sgtf=='0')],c(0.5,0.25,0.75));

quantile(c(cur$tVax4[which(cur$totVax=='4'&cur$sgtf=='1')],cur$tVax5[which(cur$totVax=='5'&cur$sgtf=='1')]),c(0.5,0.25,0.75));
quantile(c(cur$tVax4[which(cur$totVax=='4'&cur$sgtf=='0')],cur$tVax5[which(cur$totVax=='5'&cur$sgtf=='0')]),c(0.5,0.25,0.75));

table(cur$testDate,cur$sgtf)
#### sgtf
########### ba1 SGTF, delta SG
########### ba2 SG, ba1 SGTF
########### ba45 STGTF, ba2 SG
########### xbb SG, oth SGTF

table(completeDat[[1]]$sgtf,completeDat[[1]]$testWeek)

library(MASS)

set.seed(1)
parsVax = parsNat = c()
for (m in 1:10){

  cur = completeDat[[m]]

  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'

  mod = glm(as.numeric(xbb)~totVax+prior2plus+
              #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
              #as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]
  
  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))[,2:6]
  parsVax = rbind(parsVax,parsTemp)
  
  mod = glm(as.numeric(xbb)~prior2plus,#+#totVax+
              #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
              #as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]
  
  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))[,2:3]
  parsNat = rbind(parsNat,parsTemp)
  
  print(m)
}

for (j in 1:2){print(exp(quantile(parsNat[,j],c(0.5,0.025,0.975))))}
for (j in 1:5){print(exp(quantile(parsVax[,j],c(0.5,0.025,0.975))))}

table(cur$xbb)
tab = table(cur$totVax,cur$xbb); tab
for (j in 1:2){print(round(100*tab[,j]/sum(tab[,j]),1))}

types = c(cur$vaxType1,cur$vaxType2,cur$vaxType3,cur$vaxType4,cur$vaxType5)
tab = table(types[which(types!='')])
cbind(tab,round(100*tab/sum(tab),1))





set.seed(1)
parsVax = c()
for (m in 1:10){
  
  cur = completeDat[[m]]
  
  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'
  cur = cur[which(cur$numPrior!='0'),]
  
  mod = glm(as.numeric(xbb)~as.numeric(totVax)+
              as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
              as.factor(testWeek),
              data=cur,family='binomial')#[[12]][2:10,]
  
  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))[,2:6]
  parsVax = rbind(parsVax,parsTemp)
  
  print(m)
}

for (j in 1:5){print(exp(quantile(parsVax[,j],c(0.5,0.025,0.975))))}





set.seed(1)
parsImm = c()
for (m in 1:10){
  
  cur = completeDat[[m]]

  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'

  
  
  cur$imm = (cur$othImm=='1')|(cur$organ=='1')|(cur$hiv=='1')|(cur$cancer%in%c('1','carc'))|(cur$rheum=='1')
  #cur$immNum = (cur$othImm=='1')+(cur$organ=='1')+(cur$hiv=='1')+(cur$cancer%in%c('1','carc'))+(cur$rheum=='1')

  
  mod = glm(as.numeric(xbb)~imm,#+#as.factor(cancer%in%c('1','car')),#+
              #prior2plus+totVax+
              #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+bmiClass+#as.factor(chargrp)+
              #as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]
  summary(mod)
  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))[,2]
  parsImm = rbind(parsImm,parsTemp)
  
  print(m)
}
for (i in 1){print(exp(quantile(parsImm[,i],c(0.5,0.025,0.975))))}
#table(cur$sgtf)
#tab = table(imm=cur$cancer%in%c('1','carc'),xbb=cur$xbb); tab[2,]; 100*tab[2,]/table(cur$xbb)
#tab = table(imm=cur$organ,xbb=cur$xbb); tab[2,]; 100*tab[2,]/table(cur$xbb)



sum.na = function(x){return(sum(x,na.rm=T))}
for (m in 1:10){
  totPreVax = cbind(completeDat[[m]]$tInf1>=completeDat[[m]]$tVax2,completeDat[[m]]$tInf2>=completeDat[[m]]$tVax2,completeDat[[m]]$tInf3>=completeDat[[m]]$tVax2,completeDat[[m]]$tInf4>=completeDat[[m]]$tVax2)
  totPreVax = apply(totPreVax,1,sum.na)
  sel = which(completeDat[[m]]$numPrior>0&completeDat[[m]]$totVax<2)
  totPreVax[sel] = completeDat[[m]]$numPrior[sel]
  preVaxInf = totPreVax>0
  completeDat[[m]]$totPreVax = totPreVax
  completeDat[[m]]$preVaxInf = preVaxInf
  
  totPostVax = cbind(completeDat[[m]]$tInf1<completeDat[[m]]$tVax2,completeDat[[m]]$tInf2<completeDat[[m]]$tVax2,completeDat[[m]]$tInf3<completeDat[[m]]$tVax2,completeDat[[m]]$tInf4<completeDat[[m]]$tVax2)
  totPostVax = apply(totPostVax,1,sum.na)
  postVaxInf = totPostVax>0
  completeDat[[m]]$totPostVax = totPostVax
  completeDat[[m]]$postVaxInf = postVaxInf
  
  totPostInf = cbind(completeDat[[m]]$tVax1<completeDat[[m]]$tInf1,
                     completeDat[[m]]$tVax2<completeDat[[m]]$tInf1,
                     completeDat[[m]]$tVax3<completeDat[[m]]$tInf1,
                     completeDat[[m]]$tVax4<completeDat[[m]]$tInf1,
                     completeDat[[m]]$tVax5<completeDat[[m]]$tInf1)
  
  totPostInf = apply(totPostInf,1,sum.na)
  totPreInf = as.numeric(completeDat[[m]]$totVax)-totPostInf
  totPostInf[totPostInf<0] = 0; totPreInf[totPreInf<0] = 0

  completeDat[[m]]$totPostInf = totPostInf
  completeDat[[m]]$totPreInf = totPreInf
  
}



set.seed(1)
parsImm = c()
for (m in 1:10){

  cur = completeDat[[m]]
  
  cur$xbb = cur$sgtf=='0'
  cur$post2 = cur$totPostVax; cur$post2[which(cur$post2>=2)] = '2'
  cur$pre2 = cur$totPreVax; cur$pre2[which(cur$pre2>=2)] = '2'
  
  cur$post1 = cur$totPostVax; cur$post1[which(cur$post1>=2)] = '1'
  cur$pre1 = cur$totPreVax; cur$pre1[which(cur$pre1>=2)] = '1'
  
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'
  
  cur$preInf4 = cur$totPreInf; cur$preInf4[which(cur$preInf4>=4)] = '4'
  cur$postInf1 = cur$totPostInf; cur$postInf1[which(cur$postInf1>=1)] = '1'
  
  cur$vaxTiming = paste(cur$totPreInf,cur$postInf1,sep='-')

  ####
  
  #mod = glm(as.numeric(xbb)~pre1+post1+
  #            #as.factor(totVax)+
  #            #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
  #            as.factor(testWeek),
  #          data=cur,family='binomial')#[[12]][2:10,]
 
 # mod = glm(as.numeric(xbb)~totVax,#prior2plus+
 #           #  as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
 #           #  as.factor(testWeek),
 #           data=cur,family='binomial')#[[12]][2:10,]
   
  mod = glm(as.numeric(xbb)~vaxTiming+#pre2+post2+
             # as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
              as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]

  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))[,2:11]
  parsImm = rbind(parsImm,parsTemp)
  
  print(m)
}
#summary(mod)
for (j in 1:10){print(exp(quantile(parsImm[,j],c(0.5,0.025,0.975))))}


tab=table(cur$totVax,cur$xbb); tab
for (i in 1:2){print(cbind(tab[,i],100*tab[,i]/sum(tab[,i])))}


tab=table(cur$pre1,cur$xbb); tab
for (i in 1:2){print(100*tab[,i]/sum(tab[,i]))}


set.seed(1)
parsImm = c()
for (m in 1:10){

  cur = completeDat[[m]]
  
  cur$xbb = cur$sgtf=='0'
  #cur$priorOmicron = cur$priorBa12=='TRUE'
  cur$alphaEps = cur$priorAlpha=='TRUE'|cur$priorEpsilon=='TRUE'#|cur$priorWuhan=='TRUE'
  cur$preDelta = cur$priorAlpha=='TRUE'|cur$priorEpsilon=='TRUE'|cur$priorWuhan=='TRUE'

  mod = glm(as.numeric(xbb)~alphaEps,#+#priorBa45+priorBa12+priorDelta+alphaEps+
              #as.factor(totVax)+
              #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
              #as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]

  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))#[,c(2:6)]
  parsImm = rbind(parsImm,parsTemp)
  
  print(m)
}
exp(quantile(parsImm[,2],c(0.5,0.025,0.975)))

tab = table(cur=cur$xbb,past=cur$numPrior)[,2]; tab; round(100*tab/c(21870,9869),1)

x = forImp$priorDelta; mean(is.na(x)); sum(is.na(x))


set.seed(1)
parsNat = c()
for (m in 1:10){

  cur = completeDat[[m]]
  #cur = cur[which(cur$totVax%in%c('0')),]
  cur$totVax4 = cur$totVax; cur$totVax4[which(cur$totVax=='5')] = '4'

  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'

  mod = glm(as.numeric(xbb)~prior2plus*totVax4,#+
              #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
               #as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]
  
  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))#[,2:3]
  parsNat = rbind(parsNat,parsTemp)

  print(m)
}

names(coef(mod))


#for (i in c(2,3)){print(quantile(exp(parsNat[,i]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+57]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+59]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+61]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+63]),c(0.5,0.025,0.975)))}

#for (i in c(2,3)){print(quantile(exp(parsNat[,i]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+19]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+21]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+23]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+25]),c(0.5,0.025,0.975)))}

#for (i in c(2,3)){print(quantile(exp(parsNat[,i]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+6]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+8]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+10]),c(0.5,0.025,0.975)))}
#for (i in c(2,3)){print(quantile(exp(parsNat[,i]+parsNat[,i+12]),c(0.5,0.025,0.975)))}


m = 1
cur = completeDat[[m]]
cur$totVax4 = cur$totVax; cur$totVax4[which(cur$totVax=='5')] = '4'
cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'
cur$xbb = cur$sgtf=='0'
expos = paste(cur$totVax4,cur$prior2plus,sep='-')
table(expos,cur$xbb)
out = c(1213,284,22); for (i in 1:3){print(out[i]/sum(out))}

for (i in 1:2){print(quantile(exp(parsNat[,i]),c(0.5,0.025,0.975)))}
tab = table(cur$numPrior,cur$xbb); 


set.seed(1)
pars = c()
for (m in 1:10){
  
  cur = completeDat[[m]]
  cur$totVax4 = cur$totVax; cur$totVax4[which(cur$totVax=='5')] = '4'
  cur$totVaxBi = paste(cur$totVax4,cur$bivalent,sep='-')
  cur$totVax90 = paste(cur$totVax4,cur$lastvax90,sep='-')
  
  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'
  
  mod = glm(as.numeric(xbb)~totVax90,#+
              #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
              #as.factor(testWeek),
            data=cur,family='binomial')#[[12]][2:10,]
  
  parsTemp = mvrnorm(1e4,coef(mod),vcov(mod))[,2:9]
  pars = rbind(pars,parsTemp)

  print(m)
}

for (i in 1:8){print(exp(quantile(pars[,i],c(0.5,0.025,0.975))))}

tab = table(cur$totVaxBi,cur$xbb)
i = 2; tab[,i]; round(100*tab[,i]/sum(tab[,i]),1)

table(paste(cur$totVax4,cur$lastvax90,sep='-'),cur$xbb)









library(survival)

parsHosp = parsHospNat = parsHospVax = parsHospVaxBiv = parsIcu = parsVent = parsDeath = c()
for (m in 1:10){
  cur = completeDat[[m]]
  cur$ind = 1:dim(cur)[1]
  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'
  
  mod = coxph(Surv(tHosp,hosp=='TRUE')~xbb+
                #prior2plus+totVax+
                #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
                strata(testWeek),data=cur);
  coefs = mvrnorm(1e4,coef(mod),vcov(mod))
  parsHosp = c(parsHosp,coefs[,1])
  
  mod = coxph(Surv(tIcu,icu=='TRUE')~xbb+
                #prior2plus+totVax+
                #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
                strata(testWeek),data=cur);
  coefs = mvrnorm(1e4,coef(mod),vcov(mod))
  parsIcu = c(parsIcu,coefs[,1])
  
  mod = coxph(Surv(tVent,vent=='TRUE')~xbb+
                #prior2plus+totVax+
                #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
                strata(testWeek),data=cur);
  coefs = mvrnorm(1e4,coef(mod),vcov(mod))
  parsVent = c(parsVent,coefs[,1])
  
  mod = coxph(Surv(tDeath,death=='TRUE')~xbb+
                #prior2plus+totVax+
                #as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
                strata(testWeek),data=cur);
  coefs = mvrnorm(1e4,coef(mod),vcov(mod))
  parsDeath = c(parsDeath,coefs[,1])
  
  print(m)
}

sum(cur$hosp[which(cur$xbb==1)]=='TRUE'); 1e5*sum(cur$hosp[which(cur$xbb==1)]=='TRUE')/sum(cur$tHosp[which(cur$xbb==1)])
sum(cur$hosp[which(cur$xbb==0)]=='TRUE'); 1e5*sum(cur$hosp[which(cur$xbb==0)]=='TRUE')/sum(cur$tHosp[which(cur$xbb==0)])
exp(quantile(parsHosp,c(0.5,0.025,0.975)))

sum(cur$icu[which(cur$xbb==1)]=='TRUE'); 1e5*sum(cur$icu[which(cur$xbb==1)]=='TRUE')/sum(cur$tIcu[which(cur$xbb==1)])
sum(cur$icu[which(cur$xbb==0)]=='TRUE'); 1e5*sum(cur$icu[which(cur$xbb==0)]=='TRUE')/sum(cur$tIcu[which(cur$xbb==0)])
exp(quantile(parsIcu,c(0.5,0.025,0.975)))

sum(cur$vent[which(cur$xbb==1)]=='TRUE'); 1e4*sum(cur$vent[which(cur$xbb==1)]=='TRUE')/sum(cur$tVent[which(cur$xbb==1)])
sum(cur$vent[which(cur$xbb==0)]=='TRUE'); 1e4*sum(cur$vent[which(cur$xbb==0)]=='TRUE')/sum(cur$tVent[which(cur$xbb==0)])
exp(quantile(parsVent,c(0.5,0.025,0.975)))

sum(cur$death[which(cur$xbb==1)]=='TRUE'); 1e4*sum(cur$death[which(cur$xbb==1)]=='TRUE')/sum(cur$tDeath[which(cur$xbb==1)])
sum(cur$death[which(cur$xbb==0)]=='TRUE'); 1e4*sum(cur$death[which(cur$xbb==0)]=='TRUE')/sum(cur$tDeath[which(cur$xbb==0)])
exp(quantile(parsDeath,c(0.5,0.025,0.975)))





parsNat = parsVax = c()
for (m in 1:10){

  cur = completeDat[[m]]
  cur$xbb = cur$sgtf=='0'
  cur$prior2plus = cur$numPrior; cur$prior2plus[which(cur$numPrior%in%c('3','4','5'))] = '2'
  cur$prior1plus = cur$numPrior; cur$prior1plus[which(cur$numPrior%in%c('1','2','3','4','5'))] = '1'
  cur$totVax4 = cur$totVax; cur$totVax4[which(cur$totVax=='5')] = '4'
  cur$totVax4[cur$totVax4=='1'] = '0'
  
  mod = coxph(Surv(tHosp,hosp=='TRUE')~prior1plus+totVax4+
                as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
                strata(testWeek),data=cur,subset=(xbb==1));
  coefs = mvrnorm(1e4,coef(mod),vcov(mod))
  parsNat = rbind(parsNat,coefs[,1])
  
  mod = coxph(Surv(tHosp,hosp=='TRUE')~totVax4+prior2plus+
                as.factor(agegrp)+sex+race+priorOutpt+priorEd+priorInpt+smoke+as.factor(incCat)+flu+as.factor(chargrp)+bmiClass+
                strata(testWeek),data=cur,subset=(xbb==1));
  coefs = mvrnorm(1e4,coef(mod),vcov(mod))
  parsVax = rbind(parsVax,coefs[,1:4])
  
  print(m)
}
for (i in 1:2){print(quantile(exp(parsNat),c(0.5,0.025,0.975)))}
for (i in 1:3){print(quantile(exp(parsVax[,i]),c(0.5,0.025,0.975)))}

j = '4'; sum(cur$hosp[which(cur$xbb==0&cur$totVax4==j)]=='TRUE'); 1e5*sum(cur$hosp[which(cur$xbb==0&cur$totVax4==j)]=='TRUE')/sum(cur$tHosp[which(cur$xbb==0&cur$totVax4==j)])
sum(cur$hosp[which(cur$xbb==1&cur$totVax4==j)]=='TRUE'); 1e5*sum(cur$hosp[which(cur$xbb==1&cur$totVax4==j)]=='TRUE')/sum(cur$tHosp[which(cur$xbb==1&cur$totVax4==j)])

j = '1'; sum(cur$hosp[which(cur$xbb==0&cur$prior1plus==j)]=='TRUE'); 1e5*sum(cur$hosp[which(cur$xbb==0&cur$prior1plus==j)]=='TRUE')/sum(cur$tHosp[which(cur$xbb==0&cur$prior1plus==j)])
sum(cur$hosp[which(cur$xbb==1&cur$prior1plus==j)]=='TRUE'); 1e5*sum(cur$hosp[which(cur$xbb==1&cur$prior1plus==j)]=='TRUE')/sum(cur$tHosp[which(cur$xbb==1&cur$prior1plus==j)])




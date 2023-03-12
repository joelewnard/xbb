#setwd("~/Library/CloudStorage/GoogleDrive-jlewnard@berkeley.edu/My Drive/CAP KP study/xbb/data analysis")
#load('fullObsDat.Rdata')
#load('completeDat.Rdata')
#
##### all cases (broken down as hosp cases, outpt cases)
##### sgtf/no sgtf
##### hosp admissions
#
##dim(fullObsDat)
#days = unique(ymd(fullObsDat$testDate))
#days = sort(days)
#days = days[2:which(days==ymd('2023-02-23'))]
#
#dat = completeDat[[1]]
#dat$hospDate = as.character(dat$tHosp+ymd(dat$testDate))
#fullObsDat$hospDate = as.character(fullObsDat$hospDays+ymd(fullObsDat$testDate))
#
#tot = hosp = outpt = totTf = totSGTF = totNon = propSGTF = hospSGTF = hospNon = c()
#for (i in 1:length(days)){
# tot[i] = sum(fullObsDat$testDate==as.character(days[i]))
# hosp[i] = sum(fullObsDat$testDate==as.character(days[i])&fullObsDat$hospDays<=0&fullObsDat$hosp=='TRUE')
# outpt[i] = sum(fullObsDat$testDate==as.character(days[i])&fullObsDat$hospDays>0)
# totTf[i] = sum(fullObsDat$testDate==as.character(days[i])&fullObsDat$tf=='1')
# totSGTF[i] = sum(dat$sgtf[which(dat$testDate==as.character(days[i]))]=='1')
# totNon[i] = sum(dat$sgtf[which(dat$testDate==as.character(days[i]))]=='0')
# hospSGTF[i] = sum(dat$hospDate==as.character(days[i])&dat$hosp=='TRUE'&dat$sgtf=='1')
# hospNon[i] = sum(dat$hospDate==as.character(days[i])&dat$hosp=='TRUE'&dat$sgtf=='0')
# propSGTF[i] = (totSGTF/(totNon+totSGTF))[i]
#}


months = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
daylabs = paste(months[month(ymd(days))],day(ymd(days)),'')[seq(1,length(days),7)]
plot.fn = function(y1,y2,ylim,col1,col2,yby,type=1,toplab,colpoints=NA){
  plot(1,type='n',ylim=c(0,ylim),xlim=c(0,length(days)),axes=F,ann=F)
  if (type==1){
    for (i in 1:length(days)){
      polygon(x=i+c(-1,1,1,-1)*0.4,y=rep(c(0,y1[i]),each=2)+y2[i],col=rgb(col1[1],col1[2],col1[3],col1[4]),lty=0)
      polygon(x=i+c(-1,1,1,-1)*0.4,y=rep(c(0,y2[i]),each=2),col=rgb(col2[1],col2[2],col2[3],col2[4]),lty=0)
    }      
    mtext(side=2,'Cases',cex=0.6,font=1,line=2)
  } else{
    if (is.na(colpoints)==F){
      points(y1,pch=16,cex=0.5,col=rgb(0.125,0.5,0.125,0.5))
      mtext(side=2,'Proportion',cex=0.6,font=1,line=2)    
    } else{
      points(y1,pch=16,cex=0.5,col='black')
      mtext(side=2,'Cases',cex=0.6,font=1,line=2)   
    }


  }

  box(bty='l')
  axis(side=2,at=seq(0,ylim,yby),las=1,lwd.ticks=0.5,lwd=0,cex.axis=0.6)
  axis(side=1,at=0.5+seq(1,length(days),1),las=1,lwd.ticks=0.5,lwd=0,labels=NA)
  text(y=0-0.08*ylim,x=0.5+seq(1,length(days),7),daylabs,srt=45,adj=1,xpd=T,cex=0.6)
  mtext(side=1,'Test date',cex=0.6,font=1,line=2.25)
  mtext(side=3,toplab,at=-0.4*length(days),cex=0.6,line=0.125,adj=0,font=2)
}

dim(completeDat[[1]])

pdf(file='fig1.pdf',width=6.5,height=2)
par(mar=c(3.5,3,1,0.5)); par(lwd=0.5); par(mgp=c(3,0.3,0)); par(tck=-0.02)
par(mfrow=c(1,4))
plot.fn(y1=tot-totTf,y2=totTf,ylim=2500,yby=500,col1=c(1,0,0,0.25),col2=c(0,0,1,0.5),toplab='A. Total cases')
polygon(x=c(30,30,38,38),y=2500-c(0,50,50,0),col=rgb(1,0,0,0.25),lty=0)
polygon(x=c(30,30,38,38),y=2500-c(200,250,250,200),col=rgb(0,0,1,0.4),lty=0)
text(x=rep(39,2),y=c(2475,2275),c('Non-TF assay','TF assay'),cex=0.6,adj=0)

plot.fn(y1=hosp,y2=NA,ylim=80,col1=NA,col2=NA,yby=20,type=0,toplab='B. Inpatient detections')


plot.fn(y1=totSGTF,y2=totNon,ylim=1200,yby=300,col1=c(0.75,0,1,0.25),col2=c(0.125,0.5,0.125,0.5),toplab='C. Lineage frequency')
polygon(x=c(30,30,38,38)-15,y=(2400-c(0,50,50,0))/2,col=rgb(0.75,0,1,0.25),lty=0)
polygon(x=c(30,30,38,38)-15,y=(2400-c(200,250,250,200))/2,col=rgb(0.125,0.5,0.125,0.5),lty=0)
text(x=rep(39,2)-15,y=c(2375,2175)/2,c('Non-XBB/XBB.1.5','XBB/XBB.1.5'),cex=0.6,adj=0)

plot.fn(y1=1-propSGTF,y2=NA,ylim=1,col1=NA,col2=NA,yby=0.2,type=0,toplab='D. XBB/XBB.1.5 proportion',colpoints=1)
dev.off()

table(dat$sgtf)

table(dat$sgtf[which(dat$testDate==ymd('2022-12-01'))]=='0')
table(dat$sgtf[which(dat$testDate==ymd('2023-01-23'))]=='0')




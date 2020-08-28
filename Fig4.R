colMax <- function (colData) { apply(colData, MARGIN=c(2), max)}
si <- read.table("origdata/sites.txt", stringsAsFactors = F)[-15,]
si[si$e.==1173,3] <- c( "Køemelná")
rrtw <- read.table("readydata/tw150years360grains.txt", check.names = F, stringsAsFactors = F)
rr <- read.table("readydata/samples360grains.txt", check.names = F, stringsAsFactors = F)
rr <- rr[order(rr$e., rr$age),]
dr <- unique(read.table( "origdata/target_indicators.txt", check.names = F, stringsAsFactors = F)[,1])
dr[13] <- "nap"
nap <- read.table( "origdata/NAP.txt", check.names = F, stringsAsFactors = F, header = T)[,1]
target_sum <- 360

rrtw  <- data.frame(rrtw,nap=rowSums(rrtw[,colnames(rrtw) %in% nap]), stringsAsFactors = F, check.names = F)
rr    <- data.frame(rr  ,nap=rowSums(  rr[,colnames(rr)   %in% nap]), stringsAsFactors = F, check.names = F)

dr <- dr[10:13]
doBPmla <- -50
odBPsta <- 7300
mn <-0.5

head(rr)
roA <- rr[rr$e. ==1173 , c("age", dr)]
roB <- aggregate(rrtw[rrtw$e. %in% si[si$e.!=1173,1], dr], by=list(rrtw[rrtw$e. %in% si[si$e.!=1173,1],"tw"]), FUN=mean)
roBl <- aggregate(rrtw[rrtw$e. %in% si[si$e.!=1173,1], 1], by=list(rrtw[rrtw$e. %in% si[si$e.!=1173,1],"tw"]), FUN=length)


ha <- colMax(rbind(roA[,dr], roB[,dr]))/ (target_sum/100)
d=1
ha[ha==0]<- mn/2 #min(ha[ha>0])
w <- as.integer(colMax(t(as.matrix(ha)))/mn)+1
ha4 <- w*mn
ha5 <- ha4
ha5 <- ha5*10
sir <- c((max(ha5)/3)*2,ha5,max(ha5)/2,max(ha5)/8)
dre<- c( expression("Wildgras-Typ",italic("Betula"),italic("Pinus"),"NAP" ))
sede_cary <- c(0, 825, 2025, 3375, 5025)

tiff("Fig4.tiff", height = 350, width = 300, units="mm", res=300, pointsize=25, compression="lzw")
nf <- layout(t(as.matrix( seq(1,NROW(dr)+3,1))), width=sir, height=c(5))
layout.show(nf)
par(mar=c(5,0.5,16,0))
plot(1:2,1:2,  col= 'transparent',  ylab="", xlab="", xlim=c(0, 20 ), ylim=c(odBPsta, doBPmla), yaxt="n",xaxt="n", main="", bty="n")

for(d in 1:NROW(dr)){
  tsm <- (target_sum/100)
  x <- ha[d]/tsm 
  plot(x,rr[1,"age"],pch=16,cex=0.5,col="white", ylab="", xlab= "",xlim=c(0, ha[d]), ylim=c(odBPsta, doBPmla), yaxt="n", main="", bty="n", cex.lab=2, las=2, xaxt="n")
  abline(h=sede_cary, col="gray", xpd=T, lty=2, lwd=2)
  axis(3,  at=max(x/2), dre[d], xpd=TRUE, las=2, cex.axis=2, line=0, tick=F)
  
  
  axis(side=1, at = seq(0,ha4[d], mn*10/2), labels =NA ,  tcl=par("tcl")*0.4)
  axis(side=1, at = seq(0,ha4[d], mn*10), labels = c(seq(0,ha4[d], mn*10)[-NROW(seq(0,ha4[d], mn*10))],"") , cex.axis=1, las=2)#0.75)
  
  if(d==1){
    axis(2, at=seq(0,7800, 600), labels=seq(0,7800, 600), las=2, cex.axis=1.5)
  }
  
  
  polygon(x=c(0,roA[,dr[d]]/tsm ,0) , y=c(min( roA[,"age"]),roA[,"age"],max(roA[,"age"])), col="#ff3700", border = NA)
  polygon(x=c(0,roB[,dr[d]]/tsm ,0) , y=c(min( roB[,"Group.1"]),roB[,"Group.1"],max(roB[,"Group.1"])), col= "black", border = NA)
  
  lines(roA[,dr[d]]/tsm , roA[,"age"], col= "#ff3700", lwd=2)
  lines(roB[,dr[d]]/tsm , roB[,"Group.1"], col= "black", lwd=2)

  if(d==1){
    mtext("age cal yr BP", side=2, outer=F, adj = 0.45, line=6)
  }
  if(d==3){
    mtext("%", side=1,  adj = 0, line=3)
  }
  

  
}


zbt  <- roB[,"Group.1"]
krem <- roA[,"age"]
#plot(c(rep(1,NROW(zbt)), rep(2,NROW(krem))), c(zbt, krem), pch=16,cex=0.5,col=c(rep("black",NROW(zbt)), rep("#ff3700",NROW(krem))), ylab="", xlab= "",xlim=c(0, ha[d]), ylim=c(odBPsta, doBPmla), yaxt="n", main="", bty="n", cex.lab=2, las=2, xaxt="n")
par(mar=c(5,2,16,0))
plot(x=roBl$x, y=roBl$Group.1, pch=16,cex=0.5,col="white", ylab="", xlab= "",xlim=c(0, ha[d]), ylim=c(odBPsta, doBPmla), yaxt="n", main="", bty="n", cex.lab=2, las=2, xaxt="n")
segments(x0=0, y0=roBl$Group.1, roBl$x, lwd = 12, lend=2)
abline(h=sede_cary, col="gray", xpd=T, lty=2, lwd=2)
points(rep(-5,NROW(roBl[1:27,"Group.1"])),roBl[1:27,"Group.1"],pch=15, col= "#ff3700", xpd=T, cex=1.5)
axis(3,  at=0, "data coverage", xpd=TRUE, las=2, cex.axis=2, line=0, tick=F)
axis(3,  at=c(-5), c("R"), xpd=TRUE, las=1, cex.axis=1, line=-2, tick=F, col.axis =c("#ff3700"))#, "black"))
axis(3,  at=c(4),  c("Š"), xpd=TRUE, las=1, cex.axis=1, line=-2, tick=F, col.axis =c( "black"))
mm <- 1
axis(side=1, at = seq(0,12, mm*10/2), labels =NA ,  tcl=par("tcl")*0.4)
axis(side=1, at = seq(0,12, mm*10), labels = c(seq(0,12, mm*10)) , cex.axis=1, las=2)#0.75)

text(-5,5400, srt=90, "N O   D A T A" ,col="#ff3700", xpd=T, cex=1, font = 2)
axis(side=4,line=-1.3, at=sede_cary+(diff(c(sede_cary, 7300))/2), labels = rev(c("Z1", "Z2", "Z3", "Z4", "Z5")), col.ticks = "transparent", col = "transparent",las=2, cex.axis=1.75, font=1)
mtext("# sites", side=1,  adj = 0, line=3)
dev.off()

c(0, 825, 2025, 3375, 5025)+diff(c(0, 825, 2025, 3375, 5025, 7300))/2

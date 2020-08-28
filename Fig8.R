library(reshape)
colMax <- function (colData) { apply(colData, MARGIN=c(2), max)}
si <- read.table("origdata/sites.txt", stringsAsFactors = F)[-15,]
si[si$e.==1173,3] <- c( "Skláøské Valley - ref. site")
si <- si[c(11,1:10,12:nrow(si)),]
rr <- read.table("readydata/samples360grains.txt", check.names = F, stringsAsFactors = F)
rr <- rr[order(rr$e., rr$age),]
dr <- unique(read.table( "origdata/target_indicators.txt", check.names = F, stringsAsFactors = F)[,1])

f <- 1 
d <- 1
doBPmla <- -70
odBPsta <- 7200
mn <-1

filec <- si$e.


vmx <- aggregate(rr$Cerealia, by=list(rr$e.), FUN=max)
vmx <- vmx[match(filec, vmx$Group.1),]
ha <- vmx$x
ha[ha==0]<- mn/2 #min(ha[ha>0])

w <- as.integer(colMax(t(as.matrix(ha)))/mn)+1
ha4 <- w*mn
ha5 <- ha4
ha5 <- ha5*10
sir <- c(max(ha5),ha5,max(ha5)/2, max(ha5)/5)


tiff("Fig8_Cerealia.tiff", height = 350, width = 480, units="mm", res=300, pointsize=25, compression="lzw")
nf <- layout(t(as.matrix( seq(1,nrow(vmx)+2,1))), width=sir, height=c(5))
layout.show(nf)
par(mar=c(5,0.5,23,0))
plot(1:2,1:2,  col= 'transparent',  ylab="", xlab="", xlim=c(0, 20 ), ylim=c(odBPsta, doBPmla), yaxt="n",xaxt="n", main="", bty="n")
mtext("Cerealia type", side=3,  adj = 0, outer = F, cex = 2, line=20)
for(d in 1:nrow(vmx)){   
  
  y <- rr[rr$e.==vmx[d,1],"age"]
  x <- ha[d]
  plot(rep(x, NROW(y)),y,pch=16,cex=0.5,col="white", ylab="", xlab= "",xlim=c(0, ha[d]), ylim=c(odBPsta, doBPmla), yaxt="n", main="", bty="n", cex.lab=2, las=2, xaxt="n")
  abline(h=c(0, 825, 2025, 3375, 5025), col="gray", xpd=T, lty=2, lwd=2)
  axis(3,  at=max(x/2), si[d, 3], xpd=TRUE, las=2, cex.axis=2, line=-1, tick=F)
  
  
  axis(side=1, at = seq(0,ha4[d], mn), labels =NA ,  tcl=par("tcl")*0.4)
  axis(side=1, at = seq(0,ha4[d], mn*2), labels = c(seq(0,ha4[d], mn*2)[-NROW(seq(0,ha4[d], mn*2))],"") , cex.axis=1, las=2)#0.75)

  if(d==1){
    axis(2, at=seq(0,7800, 600), labels=seq(0,7800, 600), las=2, cex.axis=1.5)
  }
  
  abline(v=0)
  polygon(x=c(0,rr[rr$e.==vmx[d,1],"Cerealia"]*10 ,0) , y=c(min( y),y,max(y)), col="gray", border = NA)
  polygon(x=c(0,rr[rr$e.==vmx[d,1],"Cerealia"] ,0) , y=c(min( y),y,max(y)), col="black", border = NA)
  
  
  if(d==1){
  mtext("age cal yr BP", side=2, outer=F, adj = 0.45, line=6, cex = 1.5)
    
  }
  if(d==9){
    mtext("%", side=1,  adj = 0, line=3, cex = 1.5)
  }
}
dev.off()


vmx <- aggregate(rr$Secale, by=list(rr$e.), FUN=max)
vmx <- vmx[match(filec, vmx$Group.1),]
ha <- vmx$x
ha[ha==0]<- mn/2 #min(ha[ha>0])

w <- as.integer(colMax(t(as.matrix(ha)))/mn)+1
ha4 <- w*mn
ha5 <- ha4
ha5 <- ha5*10
sir <- c(max(ha5),ha5,max(ha5)/2, max(ha5)/5)


tiff("Fig8_Secale.tiff", height = 350, width = 480, units="mm", res=300, pointsize=25, compression="lzw")
nf <- layout(t(as.matrix( seq(1,nrow(vmx)+2,1))), width=sir, height=c(5))
layout.show(nf)
par(mar=c(5,0.5,23,0))
plot(1:2,1:2,  col= 'transparent',  ylab="", xlab="", xlim=c(0, 20 ), ylim=c(odBPsta, doBPmla), yaxt="n",xaxt="n", main="", bty="n")
mtext("Secale", side=3,  adj = 0, outer = F, cex = 2, line=20, font = 3)
for(d in 1:nrow(vmx)){   
  
  y <- rr[rr$e.==vmx[d,1],"age"]
  x <- ha[d]
  plot(rep(x, NROW(y)),y,pch=16,cex=0.5,col="white", ylab="", xlab= "",xlim=c(0, ha[d]), ylim=c(odBPsta, doBPmla), yaxt="n", main="", bty="n", cex.lab=2, las=2, xaxt="n")
  abline(h=c(0, 825, 2025, 3375, 5025), col="gray", xpd=T, lty=2, lwd=2)
  axis(3,  at=max(x/2), si[d, 3], xpd=TRUE, las=2, cex.axis=2, line=-1, tick=F)
  
  
  axis(side=1, at = seq(0,ha4[d], mn), labels =NA ,  tcl=par("tcl")*0.4)
  axis(side=1, at = seq(0,ha4[d], mn*2), labels = c(seq(0,ha4[d], mn*2)[-NROW(seq(0,ha4[d], mn*2))],"") , cex.axis=1, las=2)#0.75)
  
  if(d==1){
    axis(2, at=seq(0,7800, 600), labels=seq(0,7800, 600), las=2, cex.axis=1.5)
  }
  
  abline(v=0)
  polygon(x=c(0,rr[rr$e.==vmx[d,1],"Secale"]*10 ,0) , y=c(min( y),y,max(y)), col="gray", border = NA)
  polygon(x=c(0,rr[rr$e.==vmx[d,1],"Secale"] ,0) , y=c(min( y),y,max(y)), col="black", border = NA)
  
  
  if(d==1){
    mtext("age cal yr BP", side=2, outer=F, adj = 0.45, line=6, cex = 1.5)
    
  }
  if(d==9){
    mtext("%", side=1,  adj = 0, line=3, cex = 1.5)
  }
}
dev.off()


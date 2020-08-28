krte <- as.matrix(read.table("readydata/tetris_plot_matrix.txt"))

pps <- data.frame(matrix(nrow=ncol(krte), ncol=2, data=unlist(strsplit(colnames(krte), "_")), byrow = T), stringsAsFactors = F)
pkre <- which(pps$X2=="01")[-1] 
pulprekre <- c(pkre[1]/2,diff(pkre)/2)

dre<- c(expression(italic("Artemisia"),
                   italic("Calluna vulgaris"  ),
                   "Cerealia type",
                   "Chenopodiaceae",
                   italic("Plantago lanceolata")*"-Typ",
                   italic("Plantago media-major")*"-Typ",
                   italic("Rumex acetosa")*" type",
                   italic("Secale"),
                   italic("Urtica") 
))

tiff("Fig5.tiff", height = 350, width = 350, units="mm", res=300, pointsize=25, compression="lzw")

par(mar=c(0,3,7,2), mfrow=c(1,1), oma=c(3,2,5,1), mgp=c(3,1,0))
image(t(krte), col=rev(gray.colors(100, start = 0, end = 0.8)), axes=F, ylab="", main="")

axis(3, at= seq(0,1,1/(ncol(krte)-1))[c(1,pkre)], rep(" Skláøské Valley - reference site",9), cex.axis=.6, line=-0.7, tick=F, las=2)
axis(3,  at=(seq(0,1,1/(ncol(krte)-1))+((1/(ncol(krte)-1))/2))[c(pkre-pulprekre, ncol(krte))], dre, las=2, cex.axis=1, line=1, tick=F)
abline(v=(seq(0,1,1/(ncol(krte)-1))+((1/(ncol(krte)-1))/2))[c(1,pkre)], col="black",  lty=2, lwd=1)
abline(v=(seq(0,1,1/(ncol(krte)-1))-((1/(ncol(krte)-1))/2))[c(1,pkre)], col="black",  lty=1, xpd=T, lwd=1)

axis(2, at=rev(seq(0,1,1/49))[seq(1,50,4)], rev(rownames(krte))[seq(1,50,4)] , las=2, cex.axis=.9)
abline(h=rev(seq(0,1,1/49)-((1/49)/2))[c(6,14,23, 34)], col="black", xpd=T, lty=2, lwd=1)


par(mgp=c(3,0.1,0))

wer<-(seq(0,1,1/(ncol(krte)-1)))[2:8]
axis(1, at=wer ,labels = rep("", NROW(wer)), las=1, cex.axis=.6,  tcl=par("tcl")*0.2) 
wer2 <- (seq(0,1,1/(ncol(krte)-1)))[seq(2,8,2)]
axis(1, at=wer2 ,labels = seq(1,7,2)[1:NROW(wer2)], las=1, cex.axis=.6,  tcl=par("tcl")*0.4)  


for(i in 1:(NROW(pkre)-1))
  
{
  wer<-(seq(0,1,1/(ncol(krte)-1)))[seq(pkre[i]+1,pkre[i+1]-2,1)]
  axis(1, at=wer ,labels = rep("", NROW(wer)), las=1, cex.axis=.6,  tcl=par("tcl")*0.2) 
  wer2 <- (seq(0,1,1/(ncol(krte)-1)))[seq(pkre[i]+1,pkre[i+1]-2,2)]
  axis(1, at=wer2 ,labels = seq(1,7,2)[1:NROW(wer2)], las=1, cex.axis=.6,  tcl=par("tcl")*0.4)  
}

wer  <-(seq(0,1,1/(ncol(krte)-1)))[seq(pkre[8]+1,ncol(krte),1)]
axis(1, at=wer ,labels = rep("", NROW(wer)), las=1, cex.axis=.6,  tcl=par("tcl")*0.2) 
wer2 <-(seq(0,1,1/(ncol(krte)-1)))[(pkre[8]+1)]
axis(1, at=wer2 ,labels = seq(1,7,2)[1:NROW(wer2)], las=1, cex.axis=.6,  tcl=par("tcl")*0.4)  

par(mgp=c(5,1,0))
axis(side=4,at=rev(seq(0,1,1/49)-((1/49)/2))[c(3,10,19, 29, 41)], labels = rev(c("Z1", "Z2", "Z3", "Z4", "Z5")), col.ticks = "transparent", col = "transparent",las=2, cex.axis=1.3, font=1)

mtext("age cal yr BP", side=2, outer=T, adj = 0.35, line=0.5)
mtext("number of sites", side=1, outer=T, adj = 0.5, line=1.5)

dev.off()


rom <- as.matrix(read.table("readydata/taxa_on_site_by_tw.txt"))
nore <- as.matrix(read.table("readydata/taxa_on_site_by_tw_nodata.txt"))
ord <- read.table("readydata/labels_to_taxa_on_site_by_tw2.txt", sep="\t", header = T)
dr <- unique(read.table( "origdata/target_indicators.txt", check.names = F, stringsAsFactors = F)[,1])
si <- read.table("origdata/sites.txt", stringsAsFactors = F)
dr <- dr[1:9]

dre<- c(expression(italic("Artemisia"),
                   "Cerealia type",
                   "Chenopodiaceae",
                   italic("Plantago lanc.")*"-Typ",
                   italic("Plantago m.-m.")*"-Typ",
                   italic("Rumex acet.")*" type",
                   italic("Secale"),
                   italic("Urtica") ,
                   italic("Calluna vulgaris"  )))

si <- si[match( as.character(ord[2:14,2]),si$file),]
si[,4] <- 1:13
ord[ord$sigle %in% seq(1,13,2),"sigle"] <- NA
tiff(paste("Fig6.tiff",sep=""), height = 400, width = 700, units="mm", res=300, pointsize=25, compression="lzw")
par(mar=c(4,3,7,2), mfrow=c(1,1), oma=c(1,2,7,1))

image(t(rom), col=rev(gray.colors(100, start = 0, end = 0.8)), axes=F, ylab="", main="")
par(new=T)
image(t(nore), col="bisque", axes=F, ylab="", main="")


axis(3,  at=(seq(0,1,1/(ncol(rom)-1))+((1/(ncol(rom)-1))/2))[seq(7, 119,14)], dre, las=2, cex.axis=1.5, line=0.5, tick=F)
abline(v=rev(seq(0,1,1/(ncol(rom)-1))+((1/(ncol(rom)-1))/2))[seq(14, 112,14)+1], col="black", xpd=F, lty=2, lwd=2)

axis(2, at=rev(seq(0,1,1/49))[seq(1,50,4)], rev(rownames(rom))[seq(1,50,4)] , las=2, cex.axis=1,  lwd=2, lwd.ticks=2)
abline(h=rev(seq(0,1,1/49)-((1/49)/2))[c(6,14,23, 34)], col="black", xpd=T, lty=2, lwd=2)
par(mgp=c(3,0.6,0))
axis(3, at = seq(0,1,(1/(ncol(rom)-1))), las=2, labels =ord$sigle, cex.axis=1 ,  tcl=par("tcl")*0.4, lwd=2, lwd.ticks=2)
axis(3, at=  seq(0,1,(1/(ncol(rom)-1))*2),labels =NA,  line=0, tick=T, las=2, lwd=2, lwd.ticks=2,  tcl=par("tcl")*0.8)

par(mgp=c(5,0.5,0))
axis(side=4,at=rev(seq(0,1,1/49)-((1/49)/2))[c(3,10,19, 29, 41)], labels = rev(c("Z1", "Z2", "Z3", "Z4", "Z5")), col.ticks = "transparent", col = "transparent",las=2, cex.axis=1.3, font=1)

legend("bottom" , paste(si$V4 ,"-", si$sigle  ), ncol = 7, bty = "n", cex=1, xpd = T, inset=-0.2, x.intersp = 0.5)

mtext("age cal yr BP", side=2, outer=T, adj = 0.35, line=0.5, cex=1.5)

dev.off()



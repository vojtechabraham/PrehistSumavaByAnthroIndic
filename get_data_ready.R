### THE DATA PREPARATION CONSISTS OF FOLLOWING STEPS
### STEP 0 determination of the target pollen sum in the dataset and width of the time window 
### STEP 1 resampling of pollen counts in samples to target sum
### STEP 2 resampling of pollen counts in time windows to target sum
### STEP 3 matrix preparation of taxa on sites by time windows
### STEP 4 matrix prepaparion for tetris plot

require(reshape)
require(stringr)
colMax <- function (colData) { apply(colData, MARGIN=c(2), max)}

refid <- 1173
source('https://github.com/vojtechabraham/pollen/pollen.R', echo=F)
ag <- read.table("origdata/agde.txt")
si <- read.table("origdata/sites.txt", stringsAsFactors = F)[-15,]
n <- read.table("origdata/numbers.txt")


n <- n[!(n$newname %in% c("aaax", "NPP")),-3]                         # get rid of NPPs 
nt  <- aggregate(n$count, by=list(n$e, n$depth, n$newname), FUN=sum)  # get rid of precision of individual pollen analysts 
colnames(nt) <- colnames(n)
nt <- nt[!(nt$e==1202 & nt$depth==40),]                               
nt[nt$e==1202,1] <- 1173                                              # merge Køemelná 1 and 2
si[si$e.==refid,3] <- c( "Køemelná - reference site")

### STEP 0
        sm <- aggregate(nt$count, by=list(nt$e, nt$depth), FUN=sum) 
        sma <- merge(ag[,c(1,3,4)], sm, by=c(1,2), all.y=T)
        sma <- sma[sma$age <= 7575,  ]
        sma <- sma[order(sma$e., sma$depth),]
        di <- data.frame(dyf=NA, site=NA)
        di <- di[-1,]
        
        for (i in 1:(NROW(si))) {
          dyf <- diff(sma[sma$e.==si[i,1], "age"])
          di <- rbind(di, data.frame(dyf=dyf, site=rep(si[i,2], NROW(dyf))) )
                                                                 }
       # di[is.na(di$dyf),]
        
        
        tiff("target_sums_tw.tiff", width=1200,height = 1000, compression = "lzw", pointsize = 23)
        par(mfrow=c(2,2), mar=c(6,6,2,2), mgp=c(4,0.7,0))
        boxplot(sma$x~sma$e., names=si$sigle, las=2, ylab="grains in samples", col="gray")
        abline(h=quantile(sma$x)[2], col="red")
        #abline(h=330, col="blue")
        #barplot((sort(sma$x)), ylim=c(100,800), las=2)
        hist(sort(sma$x), breaks=80,  main = "", xlab="grains in samples", las=2, col="gray")
        text(1200,40,paste("quantile25% =", quantile(sma$x)[2]),col="red")
        #text(1200,30,paste("median =", 330),col="blue")
        abline(v=quantile(sma$x)[2], col="red")
        #abline(v=330, col="blue")
        
        boxplot(di$dyf~di$site, las=2, ylab="years between samples", col="gray")
        abline(h=quantile(di$dyf)[4], col="blue")
        hist(di$dyf, breaks=80, main = "", xlab="years between samples", las=2, col="gray")
        abline(v=quantile(di$dyf)[4], col="blue")
        text(600,100,paste("quantile75% =", quantile(di$dyf)[4]),col="blue")
        dev.off()
        
        target_sum <- 360
        cas <- cbind( cbind(seq( -75, 7425, 150), seq(75,7575,150)))
        cs <- cas[,]

### STEP 2
        sma[,5] <- paste(sma$e., sma$depth, sep="_")
        smt <- sma[sma$x >= target_sum,  ]
        no <- cast(nt, formula = e+depth~newname)
        no[is.na(no)]<- 0
        rownames(no) <- paste(no$e, no$depth, sep="_")
        
        not <- (no[smt[,5],3:ncol(no)])
        not <- as.data.frame(not)
        not <- t(not)
        #not[1:10,1:10]
    
    ### this takes some time     
        vra <- spectra_to_target_sum(not, 100, target_sum)
    ###
        vra[is.na(vra)] <- 0
        vra <- vra[rowSums(vra[,-1])!=0,]
        rownames(vra) <- vra$taxa
        
        rr <- merge(smt,t(vra[,-1]), by.x=5, by.y = 0)
        write.table(rr, "readydata/samples360grains.txt")

### STEP 2
        nr <- melt(vra)
        lowsumsamples <- merge(nt, sma[!(sma$V5 %in% smt$V5),1:2],by=1:2)        
        lowsumsamples[,5] <- paste(lowsumsamples$e, lowsumsamples$depth, sep="_")
        lowsumsamples <- lowsumsamples[,c(3,5,4)]
        colnames(lowsumsamples) <- colnames(nr)
        #nr[nr$variable %in% lowsumsamples$variable,] # 
        nr <- rbind(lowsumsamples, nr)                          # here binding of samples defcreased to target sum and samples with lower sum than target 
        
        smtw <- age_to_tw(sma, "age", cs, 1, 2)
        smtw[,7] <- paste(smtw$e., smtw$tw, sep="_") 
        smtw[,8] <- smtw$x
        smtw[smtw$x >= target_sum,  8] <- target_sum 
        smtw <- merge(smtw, aggregate(smtw$V8, by=list(smtw$V7), FUN=sum), by.x=7, by.y = 1)
        smtwa <- smtw[smtw$x.y >= target_sum,  ]
        
        nrtw <- merge(nr, smtwa[,c("V7", "V5")], by=2, all.y = T)
        nrtw <- aggregate(nrtw$value, by=list(nrtw$V7, nrtw$taxa), FUN=sum)
        now <- cast(nrtw, formula = Group.1~Group.2)
        now[is.na(now)] <- 0
        rownames(now) <- now$Group.1
        now <- now[,-1]
        now <- as.data.frame(now)
        now <- t(now)
        min(colSums(now))           # target sum should be revealed
        #now[1:10,1:10]
        
  ### this takes some time
        vrtw <- spectra_to_target_sum(now, 100, target_sum)
  ### 
        vrtw[is.na(vrtw)] <- 0
        vrtw <- vrtw[rowSums(vrtw[,-1])!=0,]
        rownames(vrtw) <- vrtw$taxa
        
        rrtw <- merge(unique(smtwa[,c(1,2,7)]),t(vrtw[,-1]), by.x=1, by.y = 0)
        write.table(rrtw, "readydata/tw150years360grains.txt")

### STEP 3
        tsm <- (target_sum/100)
        
        rrtw <- read.table("readydata/tw150years360grains.txt", check.names = F, stringsAsFactors = F)
        dr <- unique(read.table( "origdata/target_indicators.txt", check.names = F, stringsAsFactors = F)[,1])
        dr <- dr[1:9]
        
        rra <- rrtw[,c("e.","tw", dr)]
        rra <-  melt(rra, id=c("e.", "tw"))
        
        ord <- merge(si, dr)
        ord[ord$file==si[si$e.==refid,2],2] <- paste("aaa",si[si$e.==refid,2], sep = "_") 
        ord <- ord[order(ord$y, ord$file), ]
        ord[,5] <- paste(ord$y, ord$e., sep="_")
        
        rra[,4] <- rra[,4]/tsm
        rra[,5] <- paste(rra$variable, rra$e.,  sep="_")
        rpr <- cast(rra[,c(2,4,5)], tw~V5)
        rpr <- rpr[order(-rpr$tw),]
        
        rom <-  as.matrix(rpr[,ord$V5])
        colnames(rom) <- colnames(rpr)[-1]
        rownames(rom) <- rpr$tw
        nore <- as.matrix(rpr[,ord$V5])
        
        nore[!is.na(nore)] <- 0
        nore[is.na(nore)] <- 1
        nore[nore==0 ] <- NA
        #rom[rom>0] <-1
        rom[rom==0] <-NA
        rom[is.na(rom)] <- 0
        nagr <- merge(colMax(rom), ord, by.x=0, by.y=5)
        mmx <- aggregate(nagr$x, by=list(nagr$y), FUN=max)
        mxx <- mmx$x
        rom <- rom/matrix(nrow=nrow(rom), ncol=ncol(rom), data=rep(rep(mxx, each=nrow(si)),nrow(rom)), byrow = T)
        rom[rom==0] <-NA
        
        write.table(rom, "readydata/taxa_on_site_by_tw.txt")
        write.table(nore, "readydata/taxa_on_site_by_tw_nodata.txt")
        write.table(ord, "readydata/labels_to_taxa_on_site_by_tw.txt")

### STEP 4
        
        tet <- melt(rpr, id="tw")
        tet <- tet[!is.na(tet$value),]
        tet[, ncol(tet)+1] <- tet$value
        tet[tet$V5 %in% as.character(ord[ord$e. == refid,"V5"]),ncol(tet)] <-100
        tet[, "V5" ]<- substr(tet$V5,1, unlist(gregexpr( "_", tet$V5))-1)
        addit <- merge(cbind(rpr[is.na(rpr[,2]),1] ), dr)
        addit <- data.frame(addit[,1],0,addit[,2], 100)
        colnames(addit) <- colnames(tet)
        ttt <- rbind(tet, addit)
        
        
        ttt <- ttt[order(ttt$V5, ttt$tw, -ttt$V4),]
        ttt <- ttt[,c(1,3,2,4)]
        ttt[,5] <- rownames(ttt)
        ttt <- merge(ttt, mmx, by.x=2, by.y = 1)
        ttt$value <- ttt[,"value"]/ttt$x      
        ttt <- ttt[ttt$V4!=0 ,]
        
        agr <- aggregate(ttt$V5, by=list(ttt$tw, ttt$V5),length)
        
        for(i in 1:nrow(agr)){
          ttt[ttt$V5==agr[i,"Group.2"] & ttt$tw == agr[i,"Group.1"],4] <-  seq(1, agr[i,"x"], 1)  
        }
        
        ttt<-ttt[,1:4]
        ad <- aggregate(agr$x, by=list(agr$Group.2),FUN=max )
        ad <- data.frame(V5=ad$Group.1, tw=0, value=0, V4=ad$x+1 )
        ttt<- rbind(ttt,ad)
        
        ttt[,"V5"] <- paste(ttt$V5,str_pad(ttt$V4, width=2, side="left", pad="0"), sep="_")
        krte <- cast(ttt[,1:3], formula=tw~V5)
        krte <- krte[rev(order(krte$tw)),-1]
        krte[krte==0] <-NA
        rownames(krte) <- rev(unique(ttt$tw))
        
        write.table(krte, "readydata/tetris_plot_matrix.txt")



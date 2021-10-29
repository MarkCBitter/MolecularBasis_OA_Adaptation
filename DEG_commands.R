
library(edgeR)
lk1 <- read.table("LK-1/LK-1_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk2 <- read.table("LK-2/LK-2_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk3 <- read.table("LK-3/LK-3_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk4 <- read.table("LK-4/LK-4_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk5 <- read.table("LK-5/LK-5_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk6 <- read.table("LK-6/LK-6_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk7 <- read.table("LK-7/LK-7_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk8 <- read.table("LK-8/LK-8_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk9 <- read.table("LK-9/LK-9_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk10 <- read.table("LK-10/LK-10_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk11 <- read.table("LK-11/LK-11_RAW/quant.sf", header=TRUE, as.is=TRUE)
lk12 <- read.table("LK-12/LK-12_RAW/quant.sf", header=TRUE, as.is=TRUE)
counts <- data.frame( Locus=lk1$Name )
counts$LK1 <- lk1$NumReads
counts$LK2 <- lk2$NumReads
counts$LK3 <- lk3$NumReads
counts$LK4 <- lk4$NumReads
counts$LK5 <- lk5$NumReads
counts$LK6 <- lk6$NumReads
counts$LK7 <- lk7$NumReads
counts$LK8 <- lk8$NumReads
counts$LK9 <- lk9$NumReads
counts$LK10 <- lk10$NumReads
counts$LK11 <- lk11$NumReads
counts$LK12 <- lk12$NumReads
rownames( counts ) <- counts[,1]
counts <- counts[,-1]
d <- DGEList( as.matrix( counts ) )
d <- calcNormFactors( d )
cutoff <- 12
keep <- which( apply( cpm(d), 1, function(x) sum(x>=1) >= cutoff ) )
d <- d[keep,]
group <- rep( c("A", "B", "C", "D"), each=3 )
mm <- model.matrix( ~ 0 + group )
colnames( mm ) <- c("A", "B", "C", "D")
# dim went from ~150k loci to ~48k loci! good ... ?
pdf( "MDS.pdf" )
plotMDS( d, col=as.numeric(as.factor(group)) )
dev.off()
pdf( "MVtrend.pdf" )
y <- voom( d, mm, plot=TRUE )
mtext( side=3, line=0.5, text='1-Factor Model Without Intercept' )
dev.off()
fit <- lmFit( y, mm )
pdf( "dendrogram.pdf" )
library('WGCNA')
l2c <- log2( cpm( d ) + 0.00001 )
tree <- hclust( dist( t(l2c), 'euc' ), method='average' )
dcolors <- data.frame( group=labels2colors( group ) )
plotDendroAndColors( tree, colors=labels2colors(group), main='TREE' )
dev.off()
cpms <- cpm( d )
write.table( cpms, file="cpms.tsv" )

# Contrasts:
# contrast between even and variable pH profiles (at starting pH 8.1) (D-A)  # [1] correction! (D-B)
# contrast between even and variable pH profiles (at starting pH 7.4) (C-B)  # [2] correction! (C-A)
# contrast between even and variable pH profiles (sensitive time pH 8.1)     # [3] new!        (C-B)
# contrast between even and variable pH profiles (sensitive time pH 7.4)     # [4] new!        (D-A)
# contrast between starting pH values (A+D)/2-(B+C)/2                        # [5] correction! (B+D)/2 - (A+C)/2
# contrast between sensitive time pH values (A+C)/2-(B+D)/2                  # [6] correction! (A+D)/2 - (B+C)/2
# contrast between even, different pH profiles                               # [7] new!        (B-A)
# contrast between varying, different pH profiles                            # [8] new!        (C-D)
contr <- makeContrasts( D-B, C-A, C-B, D-A, (B+D)/2-(A+C)/2, (A+D)/2-(B+C)/2, B-A, C-D, levels=colnames(coef(fit)) )
tmp <- contrasts.fit( fit, contr )
tmp <- eBayes( tmp )

tmp2 <- topTable( tmp, coef=1, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
first <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( first, file="initial_pH8.1_profileEffect.txt", row.names=F, sep="\t", quote=F )
# ~5,800 DE genes
tmp2 <- topTable( tmp, coef=2, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
second <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( second, file="initial_pH7.4_profileEffect.txt", row.names=F, sep="\t", quote=F )
# ~5,400 DE genes!
tmp2 <- topTable( tmp, coef=3, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
third <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( third, file="sensitive_pH8.1_profileEffect.txt", row.names=F, sep="\t", quote=F )
# ~1,400 DE genes
tmp2 <- topTable( tmp, coef=4, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
fourth <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( fourth, file="sensitive_pH7.4_profileEffect.txt", row.names=F, sep="\t", quote=F )
# ~8,700 DE genes
tmp2 <- topTable( tmp, coef=5, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
fifth <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( fifth, file="starting_pH_Effect.txt", row.names=F, sep="\t", quote=F )
# ~2,700 DE genes
tmp2 <- topTable( tmp, coef=6, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
sixth <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( sixth, file="sensitive_pH_Effect.txt", row.names=F, sep="\t", quote=F )
# ~2,200 DE genes
tmp2 <- topTable( tmp, coef=7, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
seventh <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( seventh, file="even_pH_effect.txt", row.names=F, sep="\t", quote=F )
# only 8 DE genes
tmp2 <- topTable( tmp, coef=8, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
eighth <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( eighth, file="varying_pH_effect.txt", row.names=F, sep="\t", quote=F )
# only 39 DE genes
# requested modeling (call this ninth contrast)
pct_tot_H <- c( 17.5, 13.5, 17.1, 0.0, 0.8, 0.0, 0.0, 0.9, 0.9, 22.9, 14.0, 17.4 )
mm2 <- model.matrix( ~pct_tot_H )
colnames( mm2 ) <- c( "Intercept", "pct_tot_H" )
pdf( "MVtrend2.pdf" )
y2 <- voom( d, mm2, plot=TRUE )
dev.off()
fit2 <- lmFit( y2, mm2 )
cf2 <- contrasts.fit( fit2, coef=2 )
tmp <- eBayes( cf2 )
tmp2 <- topTable( tmp, coef=1, sort.by='P', n=Inf )
tmp2$Gene <- rownames( tmp2 )
ninth <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
write.table( ninth, file="by.pct_tot_H.txt", row.names=F, sep="\t", quote=F )

# merge all result tables with annotation into one table
first <- first[ order( first$Gene ), ]
colnames( first ) <- paste( "first.", colnames(first), sep="" )
second <- second[ order( second$Gene ), ]
colnames( second ) <- paste( "second.", colnames(second), sep="" )
third <- third[ order( third$Gene ), ]
colnames( third ) <- paste( "third.", colnames(third), sep="" )
fourth <- fourth[ order( fourth$Gene ), ]
colnames( fourth ) <- paste( "fourth.", colnames(fourth), sep="" )
fifth <- fifth[ order( fifth$Gene ), ]
colnames( fifth ) <- paste( "fifth.", colnames(fifth), sep="" )
sixth <- sixth[ order( sixth$Gene ), ]
colnames( sixth ) <- paste( "sixth.", colnames(sixth), sep="" )
seventh <- seventh[ order( seventh$Gene ), ]
colnames( seventh ) <- paste( "seventh.", colnames(seventh), sep="" )
eighth <- eighth[ order( eighth$Gene ), ]
colnames( eighth ) <- paste( "eighth.", colnames(eighth), sep="" )
ninth <- ninth[ order( ninth$Gene ), ]
colnames( ninth ) <- paste( "ninth.", colnames(ninth), sep="" )
bigTbl <- cbind( first, second[,-1], third[,-1], fourth[,-1], fifth[,-1], sixth[,-1], seventh[,-1], eighth[,-1], ninth[,-1] )
colnames( bigTbl )[1] <- "Gene"
bigBigTbl <- merge( x=annotation, y=bigTbl, by.x="geneID", by.y="Gene" )
bigBigTbl <- cbind( bigBigTbl, cpm( d )[order(rownames(d) ), ] )
write.table( bigBigTbl, file="allContrasts.tsv", row.names=F, sep="\t", quote=F )

## GO enrichment ... mainly copying Blythe's code!
library(topGO)
anno <- read.delim("../Reference/Mgal-annotation.tsv", stringsAsFactors = F, na.strings = "-")
tmp <- strsplit(anno$GO.BiologicalProcess, split = "//|;")
GO1 <- unlist(lapply(tmp, function(x)paste(x[seq(1,length(x),2)], collapse = ",")))
GO1 <- ifelse(GO1 == "NA", NA, GO1)
tmp <- strsplit(anno$GO.MolecularFunction, split = "//|;")
GO2 <- unlist(lapply(tmp, function(x)paste(x[seq(1,length(x),2)], collapse = ",")))
GO2 <- ifelse(GO2 == "NA", NA, GO2)
tmp <- strsplit(anno$Go.CellularComponent, split = "//|;")
GO3 <- unlist(lapply(tmp, function(x)paste(x[seq(1,length(x),2)], collapse = ",")))
GO3 <- ifelse(GO3 == "NA", NA, GO3)
tmp <- cbind(GO1, GO2, GO3)
GO <- apply(tmp, 1, function(x)paste(x[which(!is.na(x))], collapse = ","))
out <- cbind(anno$geneID, GO)
drop <- which(out[,2] == "")
out <- out[-drop,]
write.table(out, file = "Gene2GO.map", sep = "\t", quote = F,
            row.names = F, col.names = F)
geneID2GO <- readMappings(file = "Gene2GO.map")
# GO Biological Process enrichment testing
for ( name in c( 'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth' ) ) {
  message( name )
  comp <- eval( parse( text=name ) )  # refer to variable with name of 'first', etc.
  colnames( comp ) <- gsub( "^[a-z]*?\\.", "", colnames( comp ) )
  geneList <- comp$P.Value
  names( geneList ) <- comp$Gene
  GOdata <- new( "topGOdata",
                 ontology = "BP",
                 allGenes = geneList,
                 geneSelectionFun = function(x)x,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO )
  # Kolmogorov-Smirnov testing
  resultsKS <- runTest( GOdata, algorithm="weight01", statistic="ks" )
  tab <- GenTable( GOdata, KS=resultsKS, topNodes=length(resultsKS@score), numChar=120 )
  tab <- tab[,c(1,2,3,6)]
  names(tab)[4] <- "Raw.P.Value"
  genes.in.term <- genesInTerm(GOdata)
  genes.in.term2 <- unlist(lapply(genes.in.term, function(x)paste(x, collapse = ", ")))
  genes.in.term2 <- cbind(names(genes.in.term), genes.in.term2)
  tmp <- match(tab$GO.ID, genes.in.term2[,1])
  tab$Genes.In.Term <- genes.in.term2[tmp,2]
  outfilename <- paste( name, '.GO.BP.txt', sep='' )
  write.table( tab, outfilename, row.names=F, quote=F, sep="\t" )
}
# GO Molecular Function enrichment testing
for ( name in c( 'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth' ) ) {
  message( name )
  comp <- eval( parse( text=name ) )  # refer to variable with name of 'first', etc.
  colnames( comp ) <- gsub( "^[a-z]*?\\.", "", colnames( comp ) )
  geneList <- comp$P.Value
  names( geneList ) <- comp$Gene
  GOdata <- new( "topGOdata",
                 ontology = "MF",
                 allGenes = geneList,
                 geneSelectionFun = function(x)x,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO )
  # Kolmogorov-Smirnov testing
  resultsKS <- runTest( GOdata, algorithm="weight01", statistic="ks" )
  tab <- GenTable( GOdata, KS=resultsKS, topNodes=length(resultsKS@score), numChar=120 )
  tab <- tab[,c(1,2,3,6)]
  names(tab)[4] <- "Raw.P.Value"
  genes.in.term <- genesInTerm(GOdata)
  genes.in.term2 <- unlist(lapply(genes.in.term, function(x)paste(x, collapse = ", ")))
  genes.in.term2 <- cbind(names(genes.in.term), genes.in.term2)
  tmp <- match(tab$GO.ID, genes.in.term2[,1])
  tab$Genes.In.Term <- genes.in.term2[tmp,2]
  outfilename <- paste( name, '.GO.MF.txt', sep='' )
  write.table( tab, outfilename, row.names=F, quote=F, sep="\t" )
}
# GO Cellular Component enrichment testing
for ( name in c( 'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth' ) ) {
  message( name )  # print to screen 
  comp <- eval( parse( text=name ) )  # refer to variable with name of 'first', etc.
  colnames( comp ) <- gsub( "^[a-z]*?\\.", "", colnames( comp ) )
  geneList <- comp$P.Value
  names( geneList ) <- comp$Gene
  GOdata <- new( "topGOdata",
                 ontology = "CC",
                 allGenes = geneList,
                 geneSelectionFun = function(x)x,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO )
  # Kolmogorov-Smirnov testing
  resultsKS <- runTest( GOdata, algorithm="weight01", statistic="ks" )
  tab <- GenTable( GOdata, KS=resultsKS, topNodes=length(resultsKS@score), numChar=120 )
  tab <- tab[,c(1,2,3,6)]
  names(tab)[4] <- "Raw.P.Value"
  genes.in.term <- genesInTerm(GOdata)
  genes.in.term2 <- unlist(lapply(genes.in.term, function(x)paste(x, collapse = ", ")))
  genes.in.term2 <- cbind(names(genes.in.term), genes.in.term2)
  tmp <- match(tab$GO.ID, genes.in.term2[,1])
  tab$Genes.In.Term <- genes.in.term2[tmp,2]
  outfilename <- paste( name, '.GO.CC.txt', sep='' )
  write.table( tab, outfilename, row.names=F, quote=F, sep="\t" )
}




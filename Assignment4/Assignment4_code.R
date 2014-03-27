install.packages("ape")
install.packages("phangorn")

if(! require(affy)){
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")
  require(affy)
}


require(ape)
require(phangorn)
require(affyPLM)

setwd("~/BS32010/Assignment4")
x<-read.dna("mustelid_coI_db.fasta",format="fasta")
d<-dist.dna(x)
write.table(as.matrix(d),"distances.csv")

tr.upgma<-upgma(d)
plot(tr.upgma)
tr.upgmar<-root(tr.upgma, outgroup="AB291075")
plot(tr.upgmar);nodelabels();add.scale.bar(length=0.001)
plot(tr.upgma,"f")
tr.nj<-nj(d)
plot(tr.nj)
par(mfrow=c(1,2))
plot(tr.nj,"f")
plot(tr.nj,"p")
par(mfrow=c(1,1))
dt.upgma<-cophenetic(tr.upgma)
dmat<-as.matrix(d)
nms<-rownames(dmat)
dt.upgma<-dt.upgma[nms, nms]
dt.upgma<-as.dist(dt.upgma)
plot(dt.upgma-d,ylab="residuals", cex=0.5,main="UPGMA")
abline(h=0,lty=3)
tr.bionj<-bionj(d)
plot(tr.bionj)
tr.upgma<-as.phylo(tr.upgma)
tr.fast<-fastme.bal(d,nni=T,spr=T,tbr=T)
plot(tr.fast)
fit<-pml(tr.upgma,as.phyDat(x))
fitGTR<-pml(tr.upgma, as.phyDat(x),k=4,inv=0.2)
fitGTR<-optim.pml(fit,T)
plot(fit)
set.seed(8)
bs<-bootstrap.pml(fitGTR,bs=100,optNni=T)
plotBS(fitGTR$tr.upgma, type="f", bs) 
#Error in plot.window(...) : need finite 'xlim' values
#In addition: Warning messages:
 # 1: In min(x) : no non-missing arguments to min; returning Inf
#2: In max(x) : no non-missing arguments to max; returning -Inf
#3: In min(x) : no non-missing arguments to min; returning Inf
#4: In max(x) : no non-missing arguments to max; returning -Inf
dt.fitGTR<-cophenetic(fitGTR$tr.upgma)

hclust(d, method = "complete", members = NULL)

#this is the original code
bs<-bootstrap.pml(fit,bs=100,optNni=T)
treeBS<-plotBS(fit$tree, type="fan", bs)
treeBS<-plotBS(fit$tree, type="p", bs)

mt<-modelTest(as.phyDat(x),G=F,I=F)
dhky<-dist.dna(x,model="GTR")

install.packages("picante")
require(picante)

orig<-evol.distinct(tr.upgma,type="fair.proportion")
orig 

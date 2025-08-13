
### Distance in High Dimensions
#library(devtools)
#install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data(tissuesGeneExpression)
dim(e)
table(tissue)
# Calculate distance between rows (samples, rather than genes)
d <- dist(t(e)) # run as 'sqrt(sum((x-y)^2))' or 'sqrt(crossprod(x-y))' for each 2 rows 
class(d)
d_m <- as.matrix(d)
d_m[3,45]
x <- e['210486_at',]
y <- e['200805_at',]
d_gene <- sqrt(sum((x-y)^2))

d=dist(t(y))
as.matrix(d)[1,2]

# Reduce 2-d to 1-d
library(MASS)
library(rafalib)
library(matrixStats)
n = 100
set.seed(1)
y=t(mvrnorm(n,c(0,0),matrix(c(1,0.95,0.95,1),2,2))) 

z1 = (y[1,]+y[2,])/2   # so A should be ([1/2 1/2] [1 -1]) by row
z2 = (y[1,]-y[2,]) 
z = rbind( z1, z2) 
A <- z %*% t(y) %*% solve(y %*% t(y)) # general approach can be used to derive transformation matrices from any linear projection
thelim <- c(-3,3)
mypar(1,2)
plot(y[1,],y[2,],xlab="Twin 1 (standardized height)",ylab="Twin 2 (standardized height)",xlim=thelim,ylim=thelim)
points(y[1,1:2],y[2,1:2],col=2,pch=16)
plot(z[1,],z[2,],xlim=thelim,ylim=thelim,xlab="Average height",ylab="Differnece in height")
points(z[1,1:2],z[2,1:2],col=2,pch=16)
# to perserve Euclidean distances between dots, we need a scalar to transform A
# Transform A to Orthogonal matrix with rows as unit vectors
# Check by A-1 %*% AT = I
A_scaled <- matrix(c(1/sqrt(2),  1/sqrt(2), # row as unique vector
                     1/sqrt(2), -1/sqrt(2)), nrow = 2, byrow = TRUE)
A_scaled %*% t(A_scaled)
z_scaled <- A_scaled %*% y

plot(y[1,],y[2,],xlab="Twin 1 (standardized height)",ylab="Twin 2 (standardized height)",xlim=thelim,ylim=thelim)
points(y[1,1:2],y[2,1:2],col=2,pch=16)
plot(z_scaled[1,],z_scaled[2,],xlim=thelim,ylim=thelim,xlab="Average height",ylab="Differnece in height")
points(z[1,1:2],z[2,1:2],col=2,pch=16)



### SVD, MDS, PCA
# Mathematically
library(rafalib)
library(MASS)
n <- 100
y <- t(mvrnorm(n,c(0,0), matrix(c(1,0.95,0.95,1),2,2)))
s <- svd(y)
round(sqrt(2) * s$u , 3)

PC1 = s$d[1]*s$v[,1] # D is a n x p diagonal M
PC2 = s$d[2]*s$v[,2]
plot(PC1,PC2,xlim=c(-3,3),ylim=c(-3,3))

# Check orthonormal of U and Vcol
#sum(s$v^2)
#a <- numeric(nrow(s$v))
#for (r in seq_len(nrow(s$v))) {
#  a[r] <- (s$v[r,1])^2
#}
#sum(as.numeric(a))


# Case in seq
library(tissuesGeneExpression)
data(tissuesGeneExpression)
set.seed(1)
ind <- sample(nrow(e),500)     # sampled 500 genes
Y <- t(apply(e[ind,],1,scale)) #standardize data for illustration

s <- svd(Y)
U <- s$u
V <- s$v
D <- diag(s$d)                 #turn it into a matrix

Yhat <- U %*% D %*% t(V)
resid <- Y - Yhat
max(abs(resid))

plot(s$d) # last few are quite close to 0 (perhaps we have some replicated columns)

k <- ncol(U)-4
Yhat <- U[,1:k] %*% D[1:k,1:k] %*% t(V[,1:k])
resid <- Y - Yhat 
max(abs(resid))

plot(s$d^2/sum(s$d^2)*100, ylab="Percent variability explained")
plot(cumsum(s$d^2)/sum(s$d^2)*100,ylab="Percent variability explained",ylim=c(0,100),type="l")
# We can approximate Y with just 95
k <- 95 ##out a possible 189
Yhat <- U[,1:k] %*% D[1:k,1:k] %*% t(V[,1:k])
resid <- Y - Yhat
boxplot(resid,ylim=quantile(Y,c(0.01,0.99)),range=0)
var(as.vector(resid))/var(as.vector(Y))    # variance we lose
1-var(as.vector(resid))/var(as.vector(Y))    # variance we maintain
1-sum(s$d[1:k]^2)/sum(s$d^2)

# Highly correlated data
# Mimic a dataset with 2 cols close to each other
m <- 100
n <- 2
x <- rnorm(m)
e <- rnorm(n*m,0,0.01)
Y <- cbind(x,x)+e      # add e (#e = 2#x) to every column
cor(Y)
Y[,1]-rowMeans(Y)
Y[,2]-rowMeans(Y)

d <- svd(Y)$d
d[1]^2/sum(d^2)       # The first col can explain most variability

library(tissuesGeneExpression)
data(tissuesGeneExpression)
s = svd(e)
#  flip the sign of each column of U and respective column in V
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips
newu= sweep(s$u,2,signflips,FUN="*")   # MARGIN: 1-rows, 2 cols
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal(s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

m = rowMeans(e)
cor(as.numeric(m), s$u[,1])

# change means does not change the distance between rows
newmeans = rnorm(nrow(e)) # random values we will add to create new means
newe = e+newmeans         # we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))
# Detrend e
y = e - rowMeans(e)
s = svd(y)
resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))
z2 <- z[1:2, ]
zn_diff <- list()
for (n in seq(2,100,1)) {
  zn <- z[1:n,]
  diff <- abs(sqrt(crossprod(e[,3]-e[,45])) - sqrt(crossprod(zn[,3] - zn[,45])) )
  zn_diff[n] <- diff
}
# Distance between sample 3 and all others
distances = sqrt(apply(e[,-3]-e[,3],2,crossprod))
z2 <- z[1:2, ]
approx_distances = sqrt(apply(z2[,-3]-z2[,3],2,crossprod))
cor(distances, approx_distances, method = "spearman")

# MDS - CMD
library(tissuesGeneExpression)
library(rafalib)
data(tissuesGeneExpression)
y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)

ftissue = factor(tissue)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)

d = dist(t(e))
mds = cmdscale(d)
cor(z[2,],mds[,2])


library(GSE5859Subset)
data(GSE5859Subset)

s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)
fgroup <- factor(sampleInfo$group)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(fgroup))
legend("topleft",levels(fgroup),col=seq_along(fgroup),pch=1)

d_1 = dist(t(geneExpression))
z_1 = cmdscale(d_1)
plot(z_1[,1],z_1[,2],col=as.numeric(fgroup))

cors <- apply(z[1:10, ], 1, function(pc) cor(pc, (as.numeric(fgroup)-1)))
month = format( sampleInfo$date, "%m")
fmonth = as.numeric(factor( month)) -1
cors <- apply(z[1:10, ], 1, function(pc) cor(pc, fmonth))

u5 <- s$u[,5]
chr <- geneAnnotation$CHR
keep <- chr != "chrUn"
u5_filtered <- u5[keep]
chr_filtered <- chr[keep]
boxplot(u5_filtered ~ chr_filtered, las=2,
        main = "5th Left Singular Vector by Chromosome",
        xlab = "Chromosome", ylab = "U[,5] value")



### Clustering
library(tissuesGeneExpression)
library(rafalib)
data(tissuesGeneExpression)

d <- dist( t(e) )
mypar()
hc <- hclust(d)
hc
plot(hc, cex=0.5)
myplclust(hc, labels=tissue, lab.col=as.fumeric(tissue),cex=0.5)
abline(h=120)
hclusters <- cutree(hc, h=120)
table(true=tissue, cluster=hclusters)
hclusters <- cutree(hc, k=8)
table(true=tissue, cluster=hclusters)

# K-means_Hierarchical clustering
library(tissuesGeneExpression)
library(rafalib)
library(genefilter)
library(gplots)
library(RColorBrewer) 
data(tissuesGeneExpression)
set.seed(1)
km <- kmeans(t(e), centers=7)
names(km)
d = dist(t(e))
mds <- cmdscale(d)

mypar(1,2)
plot(mds[,1], mds[,2]) 
plot(mds[,1], mds[,2], col=km$cluster, pch=16)
table(true=tissue,cluster=km$cluster)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
rv <- rowVars(e)
idx <- order(-rv)[1:40]
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(tissue)]
head(cbind(colnames(e),cols))
heatmap.2(e[idx,], labCol=tissue,
          trace="none", 
          ColSideColors=cols, 
          col=hmcol)

# Heatmap
library(GSE5859Subset)
library(matrixStats)
library(genefilter)
library(gplots)
library(rafalib)
library(RColorBrewer)
data(GSE5859Subset)
set.seed(10)

cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]

labcol= gsub("2005-","",sampleInfo$date)  

sds =rowMads(geneExpression)
ind = order(sds,decreasing=TRUE)[1:25]

heatmap.2(geneExpression[ind,],
          col=cols,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[ind],
          labCol=labcol,
          ColSideColors=gcol,
          key=FALSE)


cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$g )

ttest = rowttests(x,g)
sds = rowSds(x)
Indexes = list(t=order(ttest$p.value)[1:50], s=order(-sds)[1:50])
ind = as.numeric(unlist(Indexes[1]))
heatmap.2(x[ind,],
          col=cols,
          trace="none",
          scale="row",
          labCol=g,
          key=FALSE)
ind = as.numeric(unlist(Indexes[2]))
heatmap.2(x[ind,],
          col=cols,
          trace="none",
          scale="row",
          labCol=g,
          key=FALSE)
# There is no relationship between g and x but with 8,793 tests some will appear significant by chance. 
# Selecting genes with the t-test gives us a deceiving result



### Classification
# Conditional Expectation
n = 10000
set.seed(1)
men = rnorm(n,176,7) 
women = rnorm(n,162,7) 
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
pr_176 <- mean(y[x == 176] == 1) 
pr <- function(a) {
  return(mean(y[x == a] == 1))
}
pr_Y <- sapply(seq(160,178), pr)
plot(seq(160,178), pr_Y, type = 'l')
abline(a=0.5,b=0)

set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]

height_lo <- loess(Y ~ X)
predict(height_lo, data.frame(X = 168))

xs = seq(160,178)
Pr =sapply(xs,function(x0) mean(Y[X==x0]))
plot(xs,Pr)
fitted=predict(height_lo,newdata=data.frame(X=xs))
lines(xs,fitted)

### kNN, cross validation
library(GSE5859Subset)
library(caret)
library(genefilter)
library(class)
library(rafalib)
RNGkind(sample.kind = "Rounding")
data(GSE5859Subset)
# define y with groups
y = factor(sampleInfo$group)
set.seed(1)
idx <-createFolds(y, 10)
idx[[3]][2]
sapply(idx, length)
sum(sapply(idx, function(i) table(y[i])))

X = t(geneExpression)
out = which(geneAnnotation$CHR%in%c("chrX","chrY"))
X = X[,-out] # to remove sex chr

m <- 8 # select top significant genes (small pvals) 
k <- 5 # #points, parameter in kNN
ind <- idx[[2]] # as test set
pvals <- rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
ind2 <- order(pvals)[1:m]
predict <- knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
sum(predict != y[ind])

# calculate for all 10 sets each as test set 
result <- sapply(idx,function(ind){
  pvals <- rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
  ind2 <- order(pvals)[1:m]
  predict <- knn(X[-ind,ind2], X[ind,ind2], y[-ind], k=k)
  sum(predict != y[ind])
})
sum(result) / length(y)

# Select best m,k
ms <- 2^c(1:11)
ks <- seq(1,9,2)
params <- expand.grid(k=ks,m=ms)

errors <- apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result <- sapply(idx,function(ind){
    pvals <- rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 <- order(pvals)[1:m]
    predict <- knn(X[-ind,ind2], X[ind,ind2], y[-ind], k=k)
    sum(predict != y[ind])
  })
  print(sum(result) / length(y))
})
params[which.min(errors),]
errors = matrix(errors,5,11)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks))

# perform p-value filtering before cross validation
X <- geneExpression
out <- which(geneAnnotation$CHR%in%c("chrX","chrY"))
X <- X[-out,] # to remove sex chr
pvals <- rowttests(X,factor(y))$p.val
X <- t(X)
errors <- apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result <- sapply(idx,function(ind){
    ind2 <- order(pvals)[1:m]
    predict <- knn(X[-ind,ind2], X[ind,ind2], y[-ind], k=k)
    sum(predict != y[ind])
  })
  print(sum(result) / length(y))
})
params[which.min(errors),]
errors = matrix(errors,5,11)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks))

# define y with month
y = factor(as.numeric(format( sampleInfo$date, "%m")=="06"))
set.seed(1)
idx <-createFolds(y, 10)
sapply(idx, length)
sum(sapply(idx, function(i) table(y[i])))

X = t(geneExpression)
out = which(geneAnnotation$CHR%in%c("chrX","chrY"))
X = X[,-out] # to remove sex chr

ms <- 2^c(1:11)
ks <- seq(1,9,2)
params <- expand.grid(k=ks,m=ms)

errors <- apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result <- sapply(idx,function(ind){
    pvals <- rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 <- order(pvals)[1:m]
    predict <- knn(X[-ind,ind2], X[ind,ind2], y[-ind], k=k)
    sum(predict != y[ind])
  })
  print(sum(result) / length(y))
})
params[which.min(errors),]
errors = matrix(errors,5,11)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks))




### Batch effect
# Admission female vs male case
library(dagdata)
data(admissions)
print( admissions )

index_m = which(admissions$Gender==1)
accepted_m= sum(admissions$Number[index_m] * admissions$Percent[index_m]/100)
applied_m = sum(admissions$Number[index_m])

index_f = which(admissions$Gender==0)
accepted_f= sum(admissions$Number[index_f] * admissions$Percent[index_f]/100)
applied_f = sum(admissions$Number[index_f])


accept_matrix <- cbind(c(accepted_f, accepted_m), 
                       c(applied_f, applied_m))
chisq.test(accept_matrix)


H_Major <- unique(as.character(admissions$Major))
H_percent <- numeric(nrow(admissions)/2)
for (r in seq_len(nrow(admissions)/2)) {
  total_apply <- admissions$Number[r] + admissions$Number[r+6]
  total_accept <- (admissions$Number[r] * admissions$Percent[r] /100+ 
                     admissions$Number[r+6] * admissions$Percent[r+6]/100)
  H_percent[r] <- total_accept / total_apply
}
H <- data.frame(cbind(H_Major, H_percent))

cor(admissions$Number[7:12], as.numeric(H$H_percent))



# Find Batch Effects
library(qvalue)
library(genefilter)
library(Biobase)
library(GSE5859)
data(GSE5859)
geneExpression = exprs(e)
sampleInfo = pData(e)

year = format(sampleInfo$date,"%y")
eth = format(sampleInfo$ethnicity)
unique(eth)
table(year, sampleInfo$ethnicity)
month.year = format(sampleInfo$date,"%m%y")
table(month.year, sampleInfo$ethnicity)

ind <- which(year%in%c("03","04") & eth=='CEU')
g <- factor(year[ind])
pvals <- rowttests(geneExpression[,ind], g)$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05)
qvals$pi0

# Modeling Batch Effects_ComBat
library(GSE5859Subset)
library(qvalue)
library(genefilter)
data(GSE5859Subset)
sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)

g <- factor(sex)
pvals <- rowttests(geneExpression, g)$p.value
qvals <- qvalue(pvals, fdr.level=0.1)
significant_qvals <- qvals$qvalues<0.1
sum(significant_qvals[geneAnnotation$CHR %in% c('chrX','chrY')]) / sum(significant_qvals)

autosomal <- !(geneAnnotation$CHR %in% c("chrX", "chrY"))
keep <- significant_qvals & autosomal
pvals_sub <- rowttests(geneExpression[keep,], month)$p.value
mean(pvals_sub<0.05)

X <- model.matrix(~sex+month)
pvals <- sapply(seq_len(nrow(geneExpression)), function(i) {
  y = geneExpression[i,]
  fit = lm(y~X-1)       # -1: Do NOT include an intercept term, β0
  summary(fit)$coefficients['Xsex', "Pr(>|t|)"]
})
qvals <- qvalue(pvals, fdr.level=0.1)
sum(qvals$qvalues < 0.1)
significant_qvals <- qvals$qvalues<0.1
sum(significant_qvals[geneAnnotation$CHR %in% c('chrX','chrY')]) / sum(significant_qvals)

month_pvals <- sapply(seq_len(nrow(geneExpression)), function(i) {
  y = geneExpression[i,]
  fit = lm(y~X-1)       # -1: Do NOT include an intercept term, β0
  summary(fit)$coefficients['Xmonth10', "Pr(>|t|)"]
})
month_qvals <- qvalue(month_pvals, fdr.level=0.1)
sum(month_qvals$qvalues < 0.1)

# Factor analysis
library(Biobase)
library(GSE5859Subset)
library(genefilter)
library(rafalib)
library(qvalue)
data(GSE5859Subset)
y = geneExpression - rowMeans(geneExpression)
cor_matrix <- cor(y)
o <- order(sampleInfo$date)
cor_ordered <- cor_matrix[o, o]
mypar(1,2)
image(cor_matrix, main = "Sample Correlation (Original Order)")
image(cor_ordered, main = "Sample Correlation (Ordered by Date)")

s <- svd(y)
pcs <- s$v[,1:2]
month = factor( format(sampleInfo$date,"%m"))
cols = as.numeric(month)[o]
mypar(2,1)
for(i in 1:2){
  plot(pcs[o,i],col=cols,xaxt="n",xlab="")    # to remove the tick labels of the plot and add the new labels with the axis function
  label = gsub("2005-","",sampleInfo$date[o])
  axis(1,1:ncol(y),label,las=2)
}

mypar(1,1)
plot(s$d^2/sum(s$d^2),ylab="% variance explained",xlab="Principal component")
abline(a=0.1,b=0)

# PC correlate to variables
cor_PC_month <- sapply(1:10, function(i) {
  fit <- s$d[i] %*% t(s$v[,i])
  cor(as.numeric(fit), as.numeric(month))
})

sex <- factor(sampleInfo$group)
cor_PC_sex <- sapply(1:10, function(i) {
  fit <- s$d[i] %*% t(s$v[,i])
  cor(as.numeric(fit), as.numeric(sex))
})

X <- model.matrix(~sex+s$v[,1:2])
pvals <- sapply(seq_len(nrow(y)), function(i) {
  y_lm = y[i,]
  fit = lm(y_lm ~ X - 1)       # -1: Do NOT include an intercept term, β0
  summary(fit)$coefficients['Xsex1', "Pr(>|t|)"]
})
qvals <- qvalue(pvals, fdr.level=0.1)
sum(qvals$qvalues < 0.1)
significant_qvals <- qvals$qvalues<0.1
sum(significant_qvals[geneAnnotation$CHR %in% c('chrX','chrY')]) / sum(significant_qvals)

# SVA
#BiocManager::install("sva")
library(sva)
library(qvalue)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])

sex = sampleInfo$group
mod = model.matrix(~sex)
svafit = sva(geneExpression,mod)
head(svafit$sv)

for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}

y <- geneExpression-rowMeans(geneExpression)
X <- model.matrix(~sex + svafit$sv)
pvals <- sapply(seq_len(nrow(y)), function(i) {
  y_lm = y[i,]
  fit = lm(y_lm ~ X - 1)       # -1: Do NOT include an intercept term, β0
  summary(fit)$coefficients['Xsex', "Pr(>|t|)"]
})
qvals <- qvalue(pvals, fdr.level=0.1)
sum(qvals$qvalues < 0.1)
significant_qvals <- qvals$qvalues < 0.1
sum(significant_qvals[geneAnnotation$CHR %in% c('chrX','chrY')]) / sum(significant_qvals)



### Quiz week4
#BiocManager::install("bladderbatch")
library(bladderbatch)
library(genefilter)
library(qvalue)
library(sva)
data(bladderdata)
edata = exprs(bladderEset)
pheno = pData(bladderEset)

expr = edata[ ,ind]
pdata = data.frame(batch = factor(pheno$batch[ind]),
                   cancer = factor(pheno$cancer[ind]))
table(pdata$batch, pdata$cancer)

# Find DEGs without removing batch effects
pvals <- rowttests(expr, factor(pdata$cancer))$p.val
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05) / length(qvals$qvalues)

# Define batch effects
X = model.matrix(~pdata$cancer + pdata$batch)

pvals <- sapply(seq_len(nrow(expr)), function(i) {
  y = expr[i,]
  fit = lm(y~X-1)       # -1: Do NOT include an intercept term, β0
  summary(fit)$coefficients['Xpdata$batch3', "Pr(>|t|)"]  # Xpdata$cancerNormal
})
qvals <- qvalue(pvals, fdr.level=0.05)
sum(qvals$significant) / length(qvals$significant)

y <- expr - rowMeans(expr)
s <- svd(y)
prop_var <- sapply(1:10, function(n){
  s$d[n]^2 / sum(s$d^2)
})
sum(prop_var>0.05)

PC1 = s$d[1]*s$v[,1] # D is a n x p diagonal M
PC2 = s$d[2]*s$v[,2]
plot(PC1,PC2,col=factor(pdata$cancer))
legend("topleft",levels(factor(pdata$cancer)),col=seq_along(factor(pdata$cancer)),pch=1)

plot(PC1,PC2,col=factor(pdata$batch))
legend("topleft",levels(factor(pdata$batch)),col=seq_along(factor(pdata$batch)),pch=1)

abs(cor(PC1, as.numeric(factor(pdata$cancer))))

# Correct batch effects
mod = model.matrix(~pdata$cancer)
svafit = sva(y,mod)
sv <- svafit$sv
head(sv)

mod0 = model.matrix(~1, data=pdata)
fpvals = f.pvalue(expr, mod, mod0)
fqvals = qvalue(fpvals)$qvalue
mean(fqvals < 0.05)

modSv = cbind(mod,sv)
mod0Sv = cbind(mod0,sv)

fpvals = f.pvalue(expr, modSv, mod0Sv)
fqvals = qvalue(fpvals)$qvalue
mean(fqvals < 0.05)

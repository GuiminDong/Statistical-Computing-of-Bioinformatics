#part1
allData<-read.csv('GSE19804_data.csv',header=T)
allData <- allData[,-1]

# select 15 genes
geneIndex <- c(31684, 21068, 30870, 36868, 12991, 19228, 25374, 40084, 21940,  6434,  1764,  3283, 43748, 22456, 20903)
allData <- allData[geneIndex,]

# divide data into three subsets
set1Index <- c(18, 47, 24, 51, 53,  3, 29, 48, 54, 58, 56, 23, 33, 27,5, 41, 11,2, 14, 40)
set2Index <- c(36, 28, 25, 37, 57, 38, 19, 20, 10, 46, 52, 59, 50, 22,  1, 12, 34, 31,  7, 49 )
set3Index <- c(55,8, 45, 30, 39, 16,  4, 17, 42, 32, 26, 35, 13, 44,9, 43, 60, 15, 21,  6)

set1 <- t(allData[,c(set1Index,set1Index+60)])
set2 <- t(allData[,c(set2Index,set1Index+60)])
set3 <- t(allData[,c(set3Index,set1Index+60)])

# Nearest neighbors classifier(k=odd integers in (3:19))
library(class)
group <-c(rep(1,20),rep(2,20))

near <- function(train,test){
  out <- c()
  for (i in 1:9){
    k <- 1+i*2
    out[i] <- sum(group!=as.numeric(knn(train,test,as.factor(group),k=k)))}
  return(out)
}

# classification error
errorA <- (near(set1,set1))
errorB <- (near(set1,set2))
errorC <- (near(set1,set3))

# choice of k: k=5, k=5, k=5
library(ggplot2)
k <- seq(3,19,by=2)
ggplot()+geom_line(aes(k,errorA,col="result A"))+geom_line(aes(k,errorB,col="result B"))+geom_line(aes(k,errorC,col="result C"))+ylab("# of error")+ggtitle("Nearest neighbor classification ")

# validation #1 A A
sum(group!=as.numeric(knn(set1,set1,as.factor(group),k=5)))
sum(group!=as.numeric(knn(set1,set2,as.factor(group),k=5)))
sum(group!=as.numeric(knn(set1,set3,as.factor(group),k=5)))

# repetition
per.near <- function(trainId,testId){
  
  per.error <- rep(NA,5000*9)
  per.k <- rep(NA,5000)
  per.vali1or2 <- rep(NA,5000)
  per.vali3 <- rep(NA,5000)
  
  set.seed(123)
  
  for (i in 1:5000){
    
    error <- rep(NA,9)
    permute <- sample(1:60)
    train <- t(allData[,c(permute[trainId],permute[trainId]+60)])
    test <- t(allData[,c(permute[testId],permute[testId]+60)])
    
    for (j in 1:9){
      k <- j*2+1     
      error[j] <- sum(group!=as.numeric(knn(train,test,as.factor(group),k)))
    }
    
    # Distribution of classification error
    per.error[((i-1)*9+1):(i*9)] <- error
    
    # Distribution of choice of k 
    per.k[i] <- order(error)[1]*2+1
    
    # Distribution of Validation
    
    # for validation 1 and 2
    per.vali1or2[i] <- sum(group!=as.numeric(knn(train,test,as.factor(group),per.k[i])))
    # for validation 3
    train2 <- t(allData[,c(permute[1:20],permute[1:20]+60)])
    test2 <- t(allData[,c(permute[41:60],permute[41:60]+60)])
    per.vali3[i] <- sum(group!=as.numeric(knn(train2,test2,as.factor(group),per.k[i])))
  }
  
  return(list(per.error,per.k,per.vali1or2, per.vali3))
  
}

resultA <- per.near(1:20,1:20)
resultB <- per.near(1:20,21:40)
resultC <- per.near(1:20,41:60)

# Distribution of classification error
k=rep(seq(3,19,by=2),5000)
par(mfrow=c(1,1),mai=(c(1.2,1.2,1.2,1.2)))
boxplot(resultA[[1]]~k,xlab="k",ylab="# of error", main="Distribution of error for result A")
boxplot(resultB[[1]]~k,xlab="k",ylab="# of error", main="Distribution of error for result B")
boxplot(resultC[[1]]~k,xlab="k",ylab="# of error", main="Distribution of error for result C")

# Distribution of choice of k 
choicek.A <- resultA[[2]]
choicek.B <- resultB[[2]]
choicek.C <- resultC[[2]]
par(mfrow=c(1,1),mai=(c(1.2,1.2,0.8,0.5)))
boxplot(c(choicek.A,choicek.B,choicek.C)~rep(c("result A","result B","result C"),each=5000),,ylab="# of k",main="Distribution of choice of k ")

# Distribution of validation
vali1 <- resultA[[3]]
vali2 <- resultB[[3]]
vali3 <- resultB[[4]]
par(mfrow=c(1,1),mai=(c(1.2,1.2,0.8,0.5)))
boxplot(c(vali1,vali2,vali3)~rep(c("validation 1","validation 2","validation 3"),each=5000),,ylab="# of error",main="Distribution of validation")

#part 2
#about data
dataset<-read.csv('GSE19804_data.csv',header=T)
data<-dataset[,2:121]
dvar<-as.matrix(apply(data,1,var))

#Hierarchical clustering 1
#lowest variance
order(dvar)[sort(rank(dvar)>54660,decreasing=T)]
new_data_low=data[c(50924, 7692, 51118, 46770, 7002, 42955, 7980, 8513, 6958,
                    22779, 50897, 6458, 9496, 858,33565),]
x<-t(new_data_low)
d1<-dist(x)
hc1=hclust(d1)
plot(hc1,labels=c(rep(0,60),rep(1,60)))
rect.hclust(hc1,k=2)



#Hierarchical clustering 2
#highest variance
order(dvar)[sort(rank(dvar)>54660)]
new_data_high=data[c(35485, 19486, 13923, 15430, 23687, 18894, 21056, 54295,
                     19282, 15686, 32955, 14160, 41833, 15173, 23132),]
y<-t(new_data_high)
d2=dist(y)
hc2=hclust(d2)
plot(hc2,labels=c(rep(0,60),rep(1,60)))
rect.hclust(hc2,k=2)

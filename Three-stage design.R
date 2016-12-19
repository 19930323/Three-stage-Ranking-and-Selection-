#Below are the coding for one situation, when theta0<-0.2, delta1<-0.05, delta2<-0.55, K<-3.

#2-stage design
theta0<-0.2
delta1<-0.05
delta2<-0.55
K<-3
a<-function(x){
  asin(sqrt(x))
}

alpha2<- function(n1,n2,y1,y2){
  part1<-matrix(rep(0,4*(n1+1)*(n1+1)),ncol=4)
  xk<-seq(from = 0, to = n1, by = 1)
  x0<-seq(from = 0, to = n1, by = 1)
  part1[,1]<-rep(xk,(n1+1))
  part1[,2]<-rep(x0,time=c(rep((n1+1),(n1+1))))
  part1[,3]<-(a(part1[,1]/n1)-a(part1[,2]/n1))
  part1[,4]<-(1-pnorm(y2*sqrt((n1+n2)/n2)-(sqrt(2)*n1/sqrt(n2))*part1[,3]))
  
  part2<-matrix(rep(0,(2+K)*(n1+1)*(n1+1)),ncol=(2+K))
  part2[,1]<-rep(xk,(n1+1))
  for (j in 1:K+1){
    part2[,j]<- (1/(j-1))*(choose(K-1,j-2))*((dbinom(part2[,1], size=n1, prob=theta0))^(j-2))* ((pbinom((part2[,1]-1), size=n1, prob=theta0))^(K+1-j))
  }
  for (i in 1:n1^2){
    part2[i,K+2]<-sum(part2[i,])-part2[i,1]
  }
  
  
  part3 <-matrix(rep(0,1*(n1+1)*(n1+1)),ncol=1)
  part3[,1] <- (part1[,3] > (y1/sqrt(2*n1)))
  
  part4 <-matrix(rep(0,1*(n1+1)*(n1+1)),ncol=1)
  part4[,1]<-(dbinom(part1[,1], size=n1, prob=theta0))*(dbinom(part1[,2], size=n1, prob=theta0))
  
  alpha1 <-cbind(part1[,4],part2[,(K+2)],part3[,1],part4[,1])
  K*sum(apply(alpha1,1,prod))
}

beta2<- function(n1,n2,y1,y2){
  part1<-matrix(rep(0,4*(n1+1)*(n1+1)),ncol=4)
  xk<-seq(from = 0, to = n1, by = 1)
  x0<-seq(from = 0, to = n1, by = 1)
  part1[,1]<-rep(xk,(n1+1))
  part1[,2]<-rep(x0,time=c(rep((n1+1),(n1+1))))
  part1[,3]<-(a(part1[,1]/n1)-a(part1[,2]/n1))
  part1[,4]<-(1-(pnorm((y2-(sqrt(2/(n1+n2))*(n1*part1[,3]+n2*(a(theta0+delta2)-a(theta0)))))/(sqrt(n2/(n1+n2))))))
  
  part2<-matrix(rep(0,(2+K)*(n1+1)*(n1+1)),ncol=(2+K))
  part2[,1]<-rep(xk,(n1+1))
  for (j in 1:K+1){
    part2[,j]<- (1/(j-1))*(choose(K-1,j-2))*((dbinom(part1[,1], size=n1, prob=theta0+delta1))^(j-2))* ((pbinom((part1[,1]-1), size=n1, prob=theta0+delta1))^(K+1-j))
  }
  for (i in 1:n1^2){
    part2[i,K+2]<-sum(part2[i,])-part2[i,1]
  }
  
  part3 <-matrix(rep(0,1*(n1+1)*(n1+1)),ncol=1)
  part3[,1] <- (part1[,3] > (y1/sqrt(2*n1)))
  
  part4 <-matrix(rep(0,1*(n1+1)*(n1+1)),ncol=1)
  part4[,1]<-(dbinom(part1[,1], size=n1, prob=theta0+delta2))*(dbinom(part1[,2], size=n1, prob=theta0))
  
  beta <-cbind(part1[,4],part2[,(K+2)],part3[,1],part4[,1])
  sum(apply(beta,1,prod))
}

#guess n1 in (35, 40), n3 in (45, 55), y2 in (0.7, 0.9), y3 in (1.9, 2.1).
#give a reasonable range of guess each time.
n1<-seq(from=35, to=40, by=2)
n2<-seq(from=45, to=55, by=2)
y1<-seq(from=0.7,to=0.9,by=0.03)
y2<- seq(from=1.9, to=2.1,by=0.02)
ln1<-length(n1)
ln2<-length(n2)
ly1<-length(y1)
ly2<-length(y2)
stage2<-matrix(rep(0,7*ln1*ln2*ly1*ly2),ncol=7)
stage2[,1]<-rep(n1,time=c(rep(ly1*ly2*ln2,ln1)))
stage2[,2]<-rep(n2,time=c(rep(ly1*ly2*ln1,ln2)))
stage2[,3]<-rep(y1,time=c(rep(ly2*ln1*ln2,ly1)))
stage2[,4]<-rep(y2,ln1*ln2*ly1)
for (i in 1:(ln1*ln2*ly1*ly2)){
  stage2[i,5]<- beta2(stage2[i,1],stage2[i,2],stage2[i,3],stage2[i,4])
  stage2[i,6]<- alpha2(stage2[i,1],stage2[i,2],stage2[i,3],stage2[i,4])
  stage2[i,7]<- (K+1)*stage2[i,1] + stage2[i,2]*(1-pnorm(stage2[i,3]))+stage2[i,2]*(1-pnorm(stage2[i,3]-(sqrt(2*stage2[i,1]))*(a(theta0+delta1)-a(theta0))))
}
data2<-stage2[((stage2[,6]>=0.047)&(stage2[,6]<=0.054)),]
data2<-data.frame(data2)
stage_2<-data2[!duplicated(data2), ]

#3-stage design
#assign p1=0.95,dalat=0.5,n1=20, psub=0.402.
p1<-0.95
theta0<-0.2
delta1<-0.05
delta2<-0.55
K<-5
n1<-20
psub<-0.402
m<-round(psub*K)
a<-function(x){
  asin(sqrt(x))
}
alpha3<- function(n2,n3,y2,y3){
  m<-round(psub*K)
  n<- n1+n2+n3
  part1<-matrix(rep(0,4*(n2+1)*(n2+1)),ncol=4)
  xk<-seq(from = 0, to = n2, by = 1)
  x0<-seq(from = 0, to = n2, by = 1)
  part1[,1]<-rep(xk,(n2+1))
  part1[,2]<-rep(x0,time=c(rep((n2+1),(n2+1))))
  part1[,3]<-(a(part1[,1]/n2)-a(part1[,2]/n2))
  part1[,4]<-(1-pnorm(y3*sqrt(n/(n1+n3))-(((sqrt(2)*n2)/sqrt(n3+n1)))*part1[,3]))
  
  part2<-matrix(rep(0,(2+m)*(n2+1)*(n2+1)),ncol=(2+m))
  part2[,1]<-rep(xk,(n2+1))
  for (j in 1:m+1){
    part2[,j]<- (1/(j-1))*(choose(m-1,j-2))*((dbinom(part2[,1], size=n2, prob=theta0))^(j-2))* ((pbinom((part2[,1]-1), size=n2, prob=theta0))^(m+1-j))
  }
  for (i in 1:n2^2){
    part2[i,m+2]<-sum(part2[i,])-part2[i,1]
  }
  
  
  part3 <-matrix(rep(0,1*(n2+1)*(n2+1)),ncol=1)
  part3[,1] <- (part1[,3] > (y2/sqrt(2*n2)))
  
  part4 <-matrix(rep(0,1*(n2+1)*(n2+1)),ncol=1)
  part4[,1]<-(dbinom(part1[,1], size=n2, prob=theta0))*(dbinom(part1[,2], size=n2, prob=theta0))
  
  alpha1<-cbind(part1[,4],part2[,(m+2)],part3[,1],part4[,1])
  K*sum(apply(alpha1,1,prod))*p1
}

beta3<-function(n2,n3,y2,y3){
  m<-round(psub*K)
  n<-n1+n2+n3
  part1<-matrix(rep(0,4*(n2+1)*(n2+1)),ncol=4)
  xk<-seq(from = 0, to = n2, by = 1)
  x0<-seq(from = 0, to = n2, by = 1)
  part1[,1]<-rep(xk,(n2+1))
  part1[,2]<-rep(x0,time=c(rep((n2+1),(n2+1))))
  part1[,3]<-(a(part1[,1]/n2)-a(part1[,2]/n2))
  part1[,4]<-(1-pnorm(y3*sqrt(n/(n1+n3))-((((sqrt(2)*n2)/sqrt(n3+n1)))*part1[,3])-sqrt(2*(n1+n3))*(a(theta0+delta2)-a(theta0))))
  
  
  part2<-matrix(rep(0,(2+m)*(n2+1)*(n2+1)),ncol=(2+m))
  part2[,1]<-rep(xk,(n2+1))
  for (j in 1:m+1){
    part2[,j]<- (1/(j-1))*(choose(m-1,j-2))*((dbinom(part2[,1], size=n2, prob=(theta0+delta1)))^(j-2))* ((pbinom((part2[,1]-1), size=n2, prob=(delta1+theta0)))^(m+1-j))
  }
  for (i in 1:n2^2){
    part2[i,m+2]<-sum(part2[i,])-part2[i,1]
  }
  
  
  part3 <-matrix(rep(0,1*(n2+1)*(n2+1)),ncol=1)
  part3[,1] <- (part1[,3] > (y2/sqrt(2*n2)))
  
  part4 <-matrix(rep(0,1*(n2+1)*(n2+1)),ncol=1)
  part4[,1]<-(dbinom(part1[,1], size=n2, prob=theta0+delta2))*(dbinom(part1[,2], size=n2, prob=theta0))
  
  beta1 <-cbind(part1[,4],part2[,(m+2)],part3[,1],part4[,1])
  sum(apply(beta1,1,prod))*p1
}

#guess n2 in (40, 50), n3 in (65, 70), y2 in (0.6, 0.75), y3 in (2, 2.2).
#give a reasonable range of guess each time.
n2<-seq(from=40, to=50, by=1)
n3<-seq(from=65, to=70, by=1)
y2<-seq(from=0.6,to=0.75,by=0.03)
y3<- seq(from=2.0, to=2.2,by=0.02)
ln2<-length(n2)
ln3<-length(n3)
ly2<-length(y2)
ly3<-length(y3)
stage3<-matrix(rep(0,7*ln2*ln3*ly2*ly3),ncol=7)
stage3[,1]<-rep(n2,time=c(rep(ly2*ly3*ln3,ln2)))
stage3[,2]<-rep(n3,time=c(rep(ly2*ly3*ln2,ln3)))
stage3[,3]<-rep(y2,time=c(rep(ly3*ln3*ln2,ly2)))
stage3[,4]<-rep(y3,ln2*ln3*ly2)
for (i in 1:(ln3*ln2*ly3*ly2)){
  stage3[i,5]<- beta3(stage3[i,1],stage3[i,2],stage3[i,3],stage3[i,4])
  stage3[i,6]<- alpha3(stage3[i,1],stage3[i,2],stage3[i,3],stage3[i,4])
  stage3[i,7]<- (K+1)*n1+(m+1)*+stage3[i,2]*(1-pnorm(stage3[i,3]))+stage3[i,2]*(1-pnorm(stage3[i,3]-(sqrt(2*stage3[i,1]))*(a(theta0+delta2)-a(theta0))))
}
data3<-stage3[((stage3[,6]>=0.049)&(stage3[,6]<=0.051)),]
data3<-data.frame(data3)
stage_3<-data2[!duplicated(data3), ]
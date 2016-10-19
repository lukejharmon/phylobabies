# pure R script to replicate full Bokma paper

# simulate trees with 19 species and root age of 46.7 mya
library(TreeSim)


# looking at Bokma's file you can deduce the following steps:
# draw net diversification rate from uniform (0.01, 0.06)
# draw extinction fraction from uniform (0, 0.9)
# back-calculate speciation and extinction rates
# make a tree
# he is using rejection sampling but I should be able to do it faster using stadler magic.

r = runif(10000, min=0.01, max=0.06)
eps = runif(10000, min=0, max=0.9)
lambda = r / (1-eps)
mu = lambda - r

simtrees<-list()
for(i in 1:10000) {
  simtrees[[i]]<-sim.bd.taxa.age(n=19, age=46.7, numbsim=1, lambda = lambda[i], mu = mu[i], mrca=T)
}

btMatrix<-matrix(nrow=10000, ncol=18)
for(i in 1:10000) {
  btMatrix[i,]<-branching.times(simtrees[[i]][[1]])
  
}


allp<-btMatrix[,-1]
allpn<-scale(log(allp))

allt = log(lambda)
mint<-min(allt)
maxt<-max(allt)
alltn<-2*(allt-min(allt))/(max(allt)-min(allt))-1

allpnpca<-prcomp(allpn)
varcomp<-allpnpca$sdev^2/sum(allpnpca$sdev^2)
allptrans<-allpnpca$x[,varcomp>0.05]
transMat<-allpnpca$rotation


# %
# %data subdivision
# p = allptrans(:,1:5000);
p = allptrans[1:5000,]


# t = alltn(1,1:5000);%only speciation rates
# something is fucked up here - I think these are extinction not speciation
t = alltn[1:5000]


# val.P = allptrans(:,5001:7000);
val.P <- allptrans[5001:7000,]

# val.T = alltn(1,5001:7000);
val.T <- alltn[5001:7000]

# cp = allptrans(:,7001:9750);
cp <- allptrans[7001:9750,1:3]

# ct = alltn(1,7001:9750);
ct <- alltn[7001:9750]

# figp = allptrans(:,9751:10001);
figp <- allptrans[9751:10000,1:3]

# figt = alltn(1,9751:10001);
figt <- alltn[9751:10000]

# %
# %training
# net = newff([minmax(p)],cells,tfunctions,'trainscg');
# net.performFcn = 'mse';
# net.trainParam.show = NaN;
# net.trainParam.lr = 0.05;
# net.trainParam.epochs = 2500;
# net.trainParam.goal = 1e-5;
# net.trainParam.max_fail = 25;
# [net,tr]=train(net,p,t,[],[],val);

dat.in <- cbind(t, p)
form.in <- as.formula(t~PC1+PC2+PC3)
seed.val <- 3
set.seed(seed.val)

# trying to mimic what Bokma did
# I don't think it matters if you use sse instead of mse for dataset of fixed size
# lr is the "learing rate"
# goal from nnet is "performance goal", while threshold from neuralnet is the stopping criteria - 
# perhaps these are the same?
# likewise I think that epochs and stepmax are the same thing
# but this function will not converge in 2500 - so I left that out.

net<-neuralnet(form.in, data=dat.in, err.fct="sse", learningrate=0.05, threshold=1e-5)

# %
# %evaluation on back-transformed data
# a = sim(net,cp);
# a = postmnmx(a,mint,maxt);%postprocess
# a = exp(a);%undo log-transform

a = compute(net, cp)$net.result
a = 0.5 * (a + 1) * (maxt - mint) + mint
a = exp(a)

# ct = postmnmx(ct,mint,maxt);%postprocess
# ct = exp(ct);
ct = 0.5 * (ct + 1) * (maxt - mint) + mint
ct = exp(ct)

# [m,b,r] = postreg(a,ct);
fm<-lm(a~ct)
m<-fm[[1]][2]
b<-fm[[1]][1]
r<-cor(ct, a)[1,1]

plot(a, ct)

# a0 = sum(ct.*a)/sum(ct.^2);%coefficient of regression through origin
a0 = sum(ct * a)/sum(ct * ct)


# save(['NetPerformance_output_ltt_' int2str(fnr) '.txt'],'m','b','r','a0','-ASCII');

res<-c(m, b, r, a0)
res



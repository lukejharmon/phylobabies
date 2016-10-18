# includes Bokma matlab code

setwd("~/Documents/machineLearning/harmon-hinchliff-learning/bokma/bokma/ANN06")

library(neuralnet)

# function y = netperformance_ltt(fnr,cells,tfunctions)
# %ltt data
# %
# %preprocessing
# load ANNvsML.txt
# ANNvsML = ANNvsML';

ANNvsML<-read.table("ANNvsML.txt")

# allp = log(ANNvsML(4:20,1:10001));
allp = log(ANNvsML[1:10001, 4:20])

# allt = log(ANNvsML(2,1:10001));
allt = log(ANNvsML[1:10001,2])

# [allpn,meanp,stdp] = prestd(allp);%normalize inputs
allpn<-scale(allp)

# [allptrans,transMat] = prepca(allpn,0.05);%principle component analysis
allpnpca<-prcomp(allpn)
varcomp<-allpnpca$sdev^2/sum(allpnpca$sdev^2)
allptrans<-allpnpca$x[,varcomp>0.05]
transMat<-allpnpca$rotation
  
# [alltn,mint,maxt] = premnmx(allt);
mint<-min(allt)
maxt<-max(allt)
alltn<-2*(allt-mint)/(maxt-mint)-1

# %
# %data subdivision
# p = allptrans(:,1:5000);
p = allptrans[1:5000,]


# t = alltn(1,1:5000);%only speciation rates
t = alltn[1:5000]


# val.P = allptrans(:,5001:7000);
val.P <- allptrans[5001:7000,]

# val.T = alltn(1,5001:7000);
val.T <- alltn[5001:7000]

# cp = allptrans(:,7001:9750);
cp <- allptrans[7001:9750,]

# ct = alltn(1,7001:9750);
ct <- alltn[7001:9750]

# figp = allptrans(:,9751:10001);
figp <- allptrans[9751:10001,]

# figt = alltn(1,9751:10001);
figt <- alltn[9751:10001]

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

# a0 = sum(ct.*a)/sum(ct.^2);%coefficient of regression through origin
a0 = sum(ct * a)/sum(ct * ct)


# save(['NetPerformance_output_ltt_' int2str(fnr) '.txt'],'m','b','r','a0','-ASCII');

res<-c(m, b, r, a0)
res



# %
# %predictions
# a = sim(net,figp);%separate data for figure and Naso
# a = postmnmx(a,mint,maxt);%postprocess
# a = exp(a);%undo log-transform

a = compute(net, figp)$net.result
a = 0.5 * (a + 1) * (maxt - mint) + mint
a = exp(a)



# a = a/a0;% (a/m)-(b/m);%reverse regression
a = a/a0

# figt = postmnmx(figt,mint,maxt);%postprocess
# figt = exp(figt);%undo log-transform
figt = 0.5 * (figt + 1) * (maxt - mint) + mint
figt = exp(figt)

plot(figt, a)

# Exc = vertcat(figt,a)';%concatenate and transpose for input to excel
# save(['NetPerformance_predictions_ltt_' int2str(fnr) '.txt'],'Exc','-ASCII');
# save (['NetPerformance_output_ltt_' int2str(fnr) '.mat'],'cells','tfunctions','net','tr');
# clear; % free some memory
# y=0;
clear all
load("data\airport.mat")
% load("urban")

%set parameters for WDCCR
gama=1e-2; 
beta=1e-2;
lamda=1e-2;
n=10; 
atomnumber=200;

%display detection map and ROC
results=hyperWDCCR(data,S,n,gama,beta,lamda,atomnumber);
figure
imagesc(results);


[pd,pf]=ROC_target(results,XY); 
pd=[0;pd;1];
pf=[0;pf;1]; 
figure
plot(pf,pd,'c');
title(f)
set(get(gca,'YLabel'),'String','Probability of detection','FontSize',10)
set(get(gca,'XLabel'),'String','False alarm rate','FontSize',10)



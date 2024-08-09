function [newT] = TargetExpasion(D,T)
%augment target pixels
% Usage
%   [newT] = expansiontarg(D,T)
% Inputs
%   D - 2d matrix of background dictionary (band*nb)
%   T - 2d matrix of target endmembers (band*ns)   
% Outputs
%   newT - target dictionary after expasion

[~,nb]=size(D);
[~,Ns] = size(T);
a=0;
newT=[];
while a<nb
    a=a+1;
    for i=0:20
        newT(:,a)=(1-i*0.01)*T(:,rem(a,Ns)+1)+(i*0.01)*D(:,a);
    end
end
end
function [results] = hyperWDCCR(X,S,n,gama,beta,lamda,atomnum)
%Perform Weighted discriminative collaborative competitive representation
%for target dection
% Usage
%   [results] = hyperWDCCR(X,S,n,gama,beta,lamda)
% Inputs
%   X - 3d matrix of HSI data (h*w*band)
%   S - 2d matrix of target endmembers (band*ns)
%   n - cluster number
%   gama,bete,lamda - parameters of WDCCR
%   atomnum - the number of background dictionary
% Outputs
%   results - vector of detector output (h*w)

atom=200;
[h,w,band]=size(X);
T=S;

maxIm=max(max(max(X)));
minIm=min(min(min(X)));
X_norm=(X-minIm)/(maxIm-minIm);
T=(T-minIm)/(maxIm-minIm);

%% removed potential target pixels
predet=TOP(X_norm,T,0.05);

%% construct the background dictionary 

D=CPSTOP(X_norm,n,predet,atomnum);

D=D';

%% target expasion
newT=TargetExpasion(D,T);
[~,Ns]=size(newT);

%% WDCCR
results=zeros(h,w);
for i=1:h
   for j=1:w
       x=X_norm(i,j,:);
       x=reshape(x,band,1); 

       %mark：the number of background dictionary
       [~,mark]=size(D);
       Wb=0;
       for nwb=1:mark
           Wb=Wb+norm(x-D(:,nwb));
       end
       Wb=Wb/mark;

       Wt=0;
       for nwt=1:Ns
           Wt=Wt+norm(x-newT(:,Ns));
       end
       Wt=Wt/Ns;

       W=blkdiag(sqrt(Wb)*eye(mark),sqrt(Wt)*eye(Ns));
       
       %get solution of WDCCR
       H=[D,newT];
       M=blkdiag(D'*D,newT'*newT);
       r=(1+gama)*(((1+2*beta)*H'*H+gama*M+lamda*(W'*W))\(H'*x));

       rb=r(1:(length(r)-Ns),1);
       rt=r((end-Ns+1):end,1);

       xx_b = D*rb;
       xx_t =newT*rt;

       results(i,j)=norm(x-xx_b)-norm(x-xx_t);

   end
end

end
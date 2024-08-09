function  [predet] = TOP(X,T,perc)
% target orthogonal purification
% Usage
%   [predet] = TOP(X,T,perc)
% Inputs
%   X - 3d matrix of HSI data (h*w*band)
%   T - 2d matrix of target endmembers (band*ns)
%   perc - the percent of removed potential target pixels
% Outputs
%   predet - the result of target orthogonal purification, binary vector(n*1) 



[h,w,band]=size(X);
X_T=orthogonalf(X,T);

X_Tback=zeros(h,w);
for i=1:h
    for j=1:w
        
        X_Tback(i,j)=dot(X_T(i,j,:),X(i,j,:));
    end
end

X_Tbackdes=sort(reshape(X_Tback,[h*w 1]),"ascend");


yita=X_Tbackdes(h*w*perc);

predet=zeros(h,w);
for i=1:h
    for j=1:w
        if X_Tback(i,j)>=yita
            predet(i,j)=1;
        else
            predet(i,j)=0;
        end
    end
end
predet=reshape(predet,[h*w 1]);
end
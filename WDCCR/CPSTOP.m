function [D] = CPSTOP(X,n,predet,atomnum)
% category-based pixel selection method with
% target orthogonal purification (CPSTOP)
% Usage
%   [D] = CPSTOP(X,n,predet,atomnum)
% Inputs
%   X - 3d matrix of HSI data (h*w*band)
%   n - the number of classes in cluster
%   predet - the result of target orthogonal purification, 2d vector
% Outputs
%   D - the background dictionary selected by CPSTOP (band*nb)


[h,w,band]=size(X);
data1=reshape(X,h*w,band);
classum=n;
opts = statset('Display','final','MaxIter',200);
idx=kmeans(data1,classum,'Options',opts);    

N=zeros(classum,1);
for i = 1:classum
    N(i) = length(find(idx == i));
end

immean=zeros(classum,band);
for i = 1:classum  
    immean(i,:) = mean(data1(find(idx == i),:));
end

D=[];

for k=1:classum
    P=round(N(k)*atomnum/(h*w));
    if(N(k)<P)
        continue
    end
    kdata=data1(find(idx==k),:);
    kpredet=predet(find(idx==k),:);
    PD=[];
    for i=1:N(k)
        PD(i)=norm(kdata(i,:)-immean(k,:));
    end
    [AscPD,index]=sort(PD,'ascend');

    
    P0=0;
    i=0;
    
    while P0<P
        i=i+1;
        if i>length(index)

            % length(index)
            break;
        end
        if kpredet(index(i))==1
            D=[D;kdata(index(i),:)];
            P0=P0+1;
        end 
    end
end

end

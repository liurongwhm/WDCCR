function [pd,pf]=ROC_target(results,XY)
% caclulate ROC of detction result
% Usage
%   [pd,pf]=ROC_target(results,XY,Outer_k)
% Inputs
%   Xresults: detection result
%   XY: the location of true target
%   Outer_k: the window size
% Outputs
%   predet - the result of target orthogonal purification, binary vector(n*1) 

[h,w,Band]=size(results);
n=h*w;
results=real(results);
N=h*w

[point,q]=size(XY);

XYRst=zeros(point,Band);
for i=1:point
   XYRst(i,:)=results( XY(i,1),XY(i,2),:); 
end
[Result]=sort(XYRst,'descend'); %


num=1:1:point;
pd=zeros(point,Band);
pf=zeros(point,Band);

for l=1:Band
    for k =1:point
            xy=zeros(h*w,2);
            mark=0;
            for i=1 : h
                for j=1 : w             
                    if ( results(i,j,l)>=Result(num(k),l) )
                        mark=mark+1;
                        xy(mark,1)=i; 
                        xy(mark,2)=j; 
                    end
                end
            end

            %
            pd_cnt=0;
            for i=1:mark
                for j=1:point
                    if(xy(i,1)==XY(j,1)&&xy(i,2)==XY(j,2))
                     pd_cnt=pd_cnt+1;
                     break;
                    end
                 end
            end
            pf_cnt=mark-pd_cnt;
            pd(k,l)=1.0*pd_cnt/point; %%Probability of detection
            pf(k,l)=pf_cnt/N;     %%False alarm rate
    end
end
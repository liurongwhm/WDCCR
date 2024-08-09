function [X_T] = orthogonalf (X,T)
%图像投影至目标正交子空间
%input: XX:待投影图像
%        T:目标向量

[h,w,band]=size(X);
XX=reshape(X,h*w,band)';


P_T = eye(band) - T * pinv(T);
X_T=P_T*XX;

%%the result after projection
X_T=reshape(X_T',[h w band]);
end
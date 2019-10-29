function [a1R] = RLS(x,namuta)
L=2;  
n = size(x,1); 
wR=zeros(L,n);
T=eye(L,L)*10;% %RLS算法下T参数的初始化,T初始值为10 
for i=(L+1):n 
    X=x(i-1:-1:(i-L));                            
    K=(T*X)/(namuta+X'*T*X);%i时刻增益值 
    e1=x(i)-wR(:,i-1)'*X; 
    wR(:,i)=wR(:,i-1)+K*e1; %i时刻权值  
    y(i)=wR(:,i)'*X; 
    e(i)=x(i)-y(i);
    T=(T-K*X'*T)/namuta; %i时刻的维纳解 
end
a1R=-wR(1,1:n); 
a2R=-wR(2,1:n); 
end

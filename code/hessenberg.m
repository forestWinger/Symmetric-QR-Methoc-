%对称矩阵的上Hessenberg化即三对角化

function [alpha,gamma,U_0]=hessenberg(A)

n=size(A,1);
U_0=eye(n,n);% 初始化三对角分解的矩阵U_0
d=zeros(n-2);

for k=1:n-2
    [v,beta]=house(A(k+1:n,k));
    u=beta*A(k+1:n,k+1:n)*v;
    w=u-(beta*u'*v/2)*v;
    A(k+1,k)=norm(A(k+1:n,k),2);
    A(k,k+1)=A(k+1,k);
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-v*w'-w*v'; %householder反射矩阵作用后的A

    d(k)=beta;
    A(k+2:n,k)=v(2:n-k);
end

v=zeros(n,1);
for k=1:n-2
    v(n+1-k:n,1)=A(n+1-k:n,n-1-k);
    v(n-k,1)=1;
    U_0=(eye(n,n)-d(n-1-k)*(v*v'))*U_0; %向后累积过程，得到U_0
end

alpha=zeros(n,1);
gamma=zeros(n-1,1);

for i=1:n
    alpha(i,1)=A(i,i);
end
for i=1:n-1
    gamma(i,1)=A(i+1,i);
end
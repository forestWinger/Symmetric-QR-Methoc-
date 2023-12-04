%观察算法1在corner case是否出现数值下溢

a=[2,4,6,8,7];
A=a'*a;
disp(A)
n=size(a,2);
u=eps;

[alpha,gamma,U_0]=hessenberg(A);

m=n; %m为检查收敛的标志
i=0; %累计迭代次数
while m>1
    [alpha(1:m),gamma(1:m-1),~]=wilkinson_QR_step(alpha(1:m),gamma(1:m-1));
    i=i+1;
    if abs(gamma(m-1))<=u*(abs(alpha(m-1))+abs(alpha(m)))
        m=m-1;
    end
    if i>=1000
        break
    end
end
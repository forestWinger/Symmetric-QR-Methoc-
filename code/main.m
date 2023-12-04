%数值算法1课程pj的主函数
%对称QR算法，用于求解对称矩阵的特征值与特征向量
%该算法应当具有O(n^3)的时间复杂度与三次收敛的特性

n_values=[5];  %更改数组得到不同结果
execution_times = zeros(size(n_values));
for i=1:size(n_values,2)
    n=n_values(i);
    history=zeros(n-1,11); %初始化收敛历史
    x=unifrnd(-100,100,n,1);
    x=sort(x);
    D=diag(x);
    U=rand(n);
    [U,R]=qr(U);
    A=U*D*U'; %初始化矩阵A，特征值为D

    tic

    [alpha,gamma,U_0]=hessenberg(A);
    Q=U_0; %初始化Q
    n=size(A,1);
    u=eps; %matlab的机器精度
    
    %算法1
%     m=n;%m为检查收敛的标志
%     while m>1
%         [alpha(1:m),gamma(1:m-1),~]=wilkinson_QR_step(alpha(1:m),gamma(1:m-1));
%         if abs(gamma(m-1))<=u*(abs(alpha(m-1))+abs(alpha(m)))
%             m=m-1;
%         end
%     end


    %改进后的算法，用于处理数值下溢
    q=0;
    j=0; %记录迭代次数
    while q<n
        % 满足条件的次对角元元素置零
        for k=1:n-1
            if abs(gamma(k))<=u*(abs(alpha(k))+abs(alpha(k+1)))
                gamma(k)=0;
            end  
        end
        
        history(1:n-1,j+1)=gamma;
        
        % 寻找最大的不可约三对角矩阵
        [p,q]=Find_Reducible(alpha,gamma);
        % QR迭代
        if q<n
            [alpha(p+1:n-q),gamma(p+1:n-q-1),G]=wilkinson_QR_step(alpha(p+1:n-q),gamma(p+1:n-q-1));
            Q(1:n,p+1:n-q)=Q(1:n,p+1:n-q)*G;
        end
        j=j+1;
    end

    execution_times(i) = toc;
end

disp(history) %显式收敛历史
    
    %绘图，用于观察时间复杂度
% figure;
% loglog(n_values, execution_times, 'o-', 'LineWidth', 2);
% xlabel('Input Size (n)');
% ylabel('Execution Time (seconds)');
% title('Algorithm Performance');
% grid on;

%误差分析
figure(1)
plot(sort(x),((sort(alpha)-sort(x)))./sort(x));
xlabel("Original Eigenvalue")
ylabel("the error of eigenvalue")
figure(2)
imagesc((Q*diag((alpha))*Q'-A)/norm(A));
colorbar
figure(3)
imagesc(Q*Q'-eye(n,n));
colorbar
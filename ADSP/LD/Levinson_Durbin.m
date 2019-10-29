%用Levinson-Durbin算法进行AR模型谱估计
%本函数实现用Levinson-Durbin算法求解Yule-Walker方程，用迭代方法求解与谱估计有关的
%参量：AR(p)阶模型回归支路的系数向量a_p，预测误差功率var_p
%===============================================================
function [a_p var_p]  = Levinson_Durbin(x,p)
%输入参量:（1）x：信号向量 （2）p：AR模型阶数
%输出参量：(1) a_p: p阶AR模型的系数 （2）var_p：p阶AR模型的预测误差功率
%首先计算自相关函数：Rxx(m) = sum(x(n)*x(n+m)/N,求和范围为从0到N-1-|m|,其中N为数据向量x的长度
%由于Matlab中数据下标从1开始，在Rxx(1)到Rxx(N)中存储x的N个自相关函数值
N = length(x);%求信号向量的长度
for ii=1:N
    Rxx(ii) = x(1:N-ii+1)*(x(ii:N))'/N;%    求自相关函数的渐近无偏估计
end
%给出迭代初始值，从1阶开始
%1阶AR模型系数a1,0=1，a1,1= -R(1)/R(0)
a(1) = 1;
a(2) = -Rxx(2)/Rxx(1);%这里由于Matlab中的下标表示a1,0 a1,1为a(1),a(2)
%用Levinson-Durbin算法进行迭代运算，迭代关系表示如下：
%var_k^2 = R(0)+sum(ak,i*R(i)),求和范围为1，2…k;
%D_k = sum(ak,i*R(k+1-i)),求和范围为从0到k，其中ak,0=1;
%gama_k+1 = D_k/var_k; var_k+1 = (1-gama_k+1^2)*var_k;
%ak+1,i = ak,i-gama_k+1*ak,k+1-i ， i=1,2…k; ak+1,k+1 = -gama_k+1;
for jj=1:p-1%迭代过程从1阶开始，直到jj+1=p阶为止
    var(jj+1) = Rxx(0+1)+a(1+1:jj+1)*Rxx(1+1:jj+1)';%预测误差功率var的计算
    D(jj+1) = a(0+1:jj+1)*(fliplr(Rxx(2:jj+1+1)))';%扩大方程中的Dk的更新
    gama(jj+1+1) = D(jj+1)/var(jj+1);%反射系数的更新
    var(jj+1+1) = (1-(gama(jj+1+1))^2)*var(jj+1);%预测误差功率var的更新
    a_temp(1) = 1;%从jj阶AR模型参数计算jj+1阶AR模型参数
    for kk=1:jj
        a_temp(kk+1) = a(kk+1)-gama(jj+1+1)*a(jj+1-kk+1);
    end
    a_temp(jj+1+1) = -gama(jj+1+1);
    a = a_temp;%将jj+1阶AR模型参数赋值给a，准备下次迭代
end
a_p = a;%迭代完成，输出p阶AR模型的回归支路系数向量
var_p = var(p+1);%输出p阶AR模型的预测误差功率
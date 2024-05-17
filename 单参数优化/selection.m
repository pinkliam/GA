%% 2.4 选择（使用轮盘赌法则）
function [newpop]=selection(pop, fitvalue) 

% 要保证fitvalue都是正数
fitvalue = fitvalue + abs(min(fitvalue)) ;

totalfit = sum(fitvalue);                       % 适应值之和
fitvalue = fitvalue./totalfit;                   % 单个个体被选择的概率
fitvalue = cumsum(fitvalue);                    % 如 fitvalue=[1 2 3 4]，则cumsum(fitvalue)=[1 3 6 10],要累加，轮盘赌法，依次看是否在转得的区域内 
[popsize,col_p] = size(pop);                     
ms = sort(rand(popsize,1));                       % 从小到大排列
fitin = 1;
newin = 1;
while newin <= popsize                           	% 选出popsize个新个体，有重复情况，和上面介绍的方法不太一样
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:) = pop(fitin,:);
        newin = newin + 1;
    else
        fitin = fitin + 1;
    end
end

end

% 选择或复制操作是决定哪些个体可以进入下一代。程序中采用赌轮盘选择法选择，这种方法较易实现。
% 根据方程 pi=fi/∑fi=fi/fsum ，选择步骤：
% 1） 在第 t 代，由（1）式计算 fsum 和 pi 
% 2） 产生 {0,1} 的随机数 rand( .)，求 s=rand( .)*fsum
% 3） 求 所有fi≥s 中最小的 k ，则第 k 个个体被选中
% 4） 进行 N 次2）、3）操作，得到 N 个个体，成为第 t=t+1 代种群
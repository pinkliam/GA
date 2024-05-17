%% 交叉
function [newpop]=crossover(pop,pc)                         % pc=0.6
[popsoze,chromlength]=size(pop);
newpop=ones(size(pop));
for i=1:2:popsoze-1                                              % 步长为2，是将相邻的两个个体进行交叉
    if(rand<pc)
        cpoint=round(rand*chromlength);
        newpop(i,:)=[pop(i,1:cpoint),pop(i+1,cpoint+1:chromlength)];
        newpop(i+1,:)=[pop(i+1,1:cpoint),pop(i,cpoint+1:chromlength)];
    else
        newpop(i,:)=pop(i,:);
        newpop(i+1,:)=pop(i+1,:);
    end
end

end

% 交叉(crossover)，群体中的每个个体之间都以一定概率 pc 交叉，即两个个体从各自字符串的某一位置
% （一般是随机确定）开始互相交换，这类似生物进化过程中的基因分裂与重组。例如，假设2个父代个体x1，x2为：
% x1=0100110
% x2=1010001
% 从每个个体的第3位开始交叉，交又后得到2个新的子代个体y1，y2分别为：
% y1＝0100001
% y2＝1010110
% 这样2个子代个体就分别具有了2个父代个体的某些特征。
% 利用交又我们有可能由父代个体在子代组合成具有更高适合度的个体。
% 交叉是遗传算法区别于其它传统优化方法的主要特点之一。


%% 交叉也可以设置成每个变量交叉

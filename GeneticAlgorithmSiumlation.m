%% 遗传算法仿真
% 功能： 
% 一元函数求极值
% Created by YanBin Wang
% Date: 2024.5.16

function GeneticAlgorithmSiumlation

clear;
clc;                        % 清屏
close all
x_max   = 10;
x_min   = 0;
x_range = [x_min, x_max];
[num_var,temp] = size(x_range);
popsize = 100;           	% 群体大小
chromlength = 20;          	% 字符串长度(个体长度)，需要根据问题求解的精度、变量的取值范围综合判定
pc = 0.7;               	% 交叉概率，只有在随机数小于pc时，才会产生交叉 一般取 60~100%
pm = 0.05;                	% 变异概率，一般 0.1~10%
iter_num    = 200;          % 遗传代数

pop = initpop(popsize, chromlength, num_var);             % 随机产生二进制编码的初始群体
for i = 1:1:iter_num              
    % 计算目标函数
    objvalue    = calobjvalue(pop, x_range, chromlength); 
    
    % 计算群体中每个个体的适应度
    fitvalue    = calfitvalue(objvalue);             	
    
    % 选择，复制
    newpop      = selection(pop,fitvalue);            	
    
    % 交叉
    newpop1     = crossover(newpop,pc);                
    
    % 变异
    newpop2     = mutation(newpop1,pm);                 
    
    % 计算目标函数
    objvalue    = calobjvalue(newpop2, x_range, chromlength); 
    
    % 计算群体中每个个体的适应度
    fitvalue    = calfitvalue(objvalue);            	
    
    % 平均适应度
    fitness_ave(i) = calavefitness(fitvalue);           
    
    % 求出群体中适应值最大的个体及其适应值
    [bestindividual,bestfit]=best(newpop2,fitvalue); 	
    
    % 返回的 y 是自适应度值，而非函数值
    y(i)    = bestfit;                                              
    
    % 将自变量解码成十进制
    x(i,:)    = decodechrom(bestindividual, x_range, chromlength);	
    
    % 更新种群
    pop     = newpop2;
end

% 作图
figure(1)
plot(1:iter_num, fitness_ave, 'r', 1:iter_num, y, 'b')
grid on
legend('平均适应度', '最优适应度')


figure(2)
fplot('x+10*sin(5*x)+7*cos(4*x)',x_range);
hold on
plot(x,y,'r*');                                       

%输出结果
[z ,index]=max(y);             %计算最大值及其位置
ymax = z;
disp(['最优染色体为', num2str(x(index,:))])
disp(['最优适应度为', num2str(ymax)])

end


%% 2.1初始化(编码)
% initpop.m函数的功能是实现群体的初始化，popsize表示群体的大小，chromlength表示染色体的长度(二值数的长度)，
% rand随机产生每个单元为 {0,1} 行数为popsize，列数为chromlength的矩阵，
% round对矩阵的每个单元进行圆整四舍五入。
% 长度大小取决于变量的二进制编码的长度*变量个数。
% 初始化
function pop=initpop(popsize,chromlength, num_var) 
if num_var == 1                 % 只有一个变量的情况
    pop = round(rand(popsize,chromlength)); 
else                            % 多个变量的情况

end

end

%% 2.2 计算目标函数值
% 2.2.1 calobjvalue.m函数的功能是实现目标函数的计算
% 实现目标函数的计算，将 二值域 中的数转化为 变量域的数
% 参数pop 表示种群数
% x_range 变量的取值范围
% chromlength  二进制编码的长度
function [objvalue]=calobjvalue(pop, x_range, chromlength)
x = decodechrom(pop, x_range, chromlength);   %将pop每行转化成十进制数
objvalue = x+10*sin(5*x)+7*cos(4*x);   %计算目标函数值

end


% 2.2.2 将二进制编码转化为十进制数(2)
% decodechrom.m函数的功能是将染色体(或二进制编码)转换为十进制，
% 参数pop 表示种群数
% x_range 变量的取值范围
% chromlength  二进制编码的长度
% (对于多个变量而言，如有两个变量，采用2*chromlength位表示，每个变量chromlength位
% 则第一个变量从1开始，另一个变量从chromlength+1开始)
function pop2=decodechrom(pop,x_range,chromlength)          %1  10

pop1 = pop;
x_min = x_range(1);
x_max = x_range(2);
pop2 = x_min + (x_max - x_min) * decodebinary(pop1)/(2^chromlength - 1);  %将pop每行转换成十进制值
end

% 将二进制数转化为十进制数clc
function pop2=decodebinary(pop)
[px,py]=size(pop);                  % 求pop行和列数
for i=1:py
    pop1(:,i)=2.^(py-i).*pop(:,i);  % 产生 [2^n 2^(n-1) ... 1] 的行向量，然后求和
end
pop2=sum(pop1,2);                   % 求pop1的每行之和，得到对应的十进制数

end

%% 2.3 计算个体的适应值
function fitvalue=calfitvalue(objvalue)
[px,py]=size(objvalue);                   %目标值有正有负
for i=1:px
    if objvalue(i)>0
        temp=objvalue(i);
    else
        temp=0.0;
    end
    fitvalue(i)=temp;
end
fitvalue=fitvalue';

end

%% 求平均适应度
function fitness_ave = calavefitness(fitness)
[N ,~] = size(fitness);
fitness_ave = sum(fitness)/N;
end

%% 2.4 选择
function [newpop]=selection(pop, fitvalue) 
totalfit = sum(fitvalue);                       % 求适应值之和
fitvalue = fitvalue/totalfit;                   % 单个个体被选择的概率
fitvalue = cumsum(fitvalue);                    % 如 fitvalue=[1 2 3 4]，则cumsum(fitvalue)=[1 3 6 10],要累加，轮盘赌法，依次看是否在转得的区域内 
[row_p,col_p] = size(pop);                      % row_p*col_p
ms = sort(rand(row_p,1));                       % 从小到大排列
fitin = 1;
newin = 1;
while newin <= row_p                           	% 选出col_p个新个体，有重复情况，和上面介绍的方法不太一样
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:) = pop(fitin,:);
        newin = newin + 1;
    else
        fitin = fitin + 1;
    end
end

end


%% 2.5 交叉
function [newpop]=crossover(pop,pc)                         % pc=0.6
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:2:px-1                                              % 步长为2，是将相邻的两个个体进行交叉
    if(rand<pc)
        cpoint=round(rand*py);
        newpop(i,:)=[pop(i,1:cpoint),pop(i+1,cpoint+1:py)];
        newpop(i+1,:)=[pop(i+1,1:cpoint),pop(i,cpoint+1:py)];
    else
        newpop(i,:)=pop(i,:);
        newpop(i+1,:)=pop(i+1,:);
    end
end

end

%% 2.6 变异
function [newpop]=mutation(pop,pm)
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:px
    if(rand<pm)
        mpoint=round(rand*py);              % 产生的变异点在1-10之间
        if mpoint<=0
            mpoint=1;                   	% 变异位置
        end
        newpop(i,:)=pop(i,:);
        if any(newpop(i,mpoint))==0
            newpop(i,mpoint)=1;
        else
            newpop(i,mpoint)=0;
        end
    else
        newpop(i,:)=pop(i,:);
    end
end

end


%% 2.7 求出群体中最大得适应值及其个体
% 遗传算法子程序
% 求出第 t 代群体中适应值最大的值
function [bestindividual,bestfit]=best(pop,fitvalue)
[px,py]=size(pop);
bestindividual=pop(1,:);
bestfit=fitvalue(1);
for i=2:px
    if fitvalue(i)>bestfit
        bestindividual=pop(i,:);
        bestfit=fitvalue(i);
    end
end

end



%% 遗传算法仿真
% Function：对多元函数求极大值(多参数)(>=1均可)
% Created by YanBin Wang
% Date: 2024.05.17

%% question
% 最优的参数为什么会出现三维显示出界情况？ 个人感觉这是显示问题，实际是没有问题
%%


warning off; % 关闭所有警告
clear;
clc;                    
close all

% 函数参数设置
max_region = 8;
min_region = 0;
lu = [ones(1,2).*min_region; ones(1,2).*max_region];
[~, num_var] = size(lu);

f = @(x, y) abs( sin(pi*(x - 3))./(pi*(x - 3)) ).*abs( sin(pi*(y - 3))./(pi*(y - 3)) );   % 目标函数

% 种群参数设置
popsize = 100;           	% 种群大小
chromlength = 20;          	% 二进制编码长度(单个变量的长度)，需要根据问题求解的精度、变量的取值范围综合判定
pc = 0.8;               	% 交叉概率，只有在随机数小于pc时，才会产生交叉 一般取 60~100%
pm = 0.07;                	% 变异概率，一般 0.01~0.1
iter_max = 200;              % 迭代次数（遗传次数）

fitness_ave = -inf(iter_max, 1); % 记录每次迭代的平均适应度
bestfit = -inf;
bestfits = -inf(iter_max, 1); % 记录最优变化

value = -inf(iter_max, num_var);
fitness = -inf(iter_max,1);

% 初始化种群，个体长度=变量的二进制编码的长度*变量个数            
pop = round(rand(popsize, chromlength*num_var)); 

% 迭代开始
for i = 1:1:iter_max             
    % 种群十进制形式
    pop_decimal = BinaryToDecimal(pop, lu, num_var);
    % 计算目标函数
    objvalue = calobjvalue(pop_decimal);
    % 计算群体中每个个体的适应度
    fitvalue = calfitvalue(objvalue);    
    % 选择，复制
    newpop_selection = selection(pop, fitvalue);            	
    % 交叉
    newpop_crossover = crossover(newpop_selection, pc);                
    % 变异
    newpop_mutation = mutation(newpop_crossover, pm);                 
    % 种群十进制形式
    pop_decimal = BinaryToDecimal(newpop_mutation, lu, num_var);
    % 计算目标函数
    objvalue = calobjvalue(pop_decimal);
    % 计算群体中每个个体的适应度
    fitvalue = calfitvalue(objvalue);  
    % 计算平均适应度
    fitness_ave(i) = calavefitness(fitvalue);  
    % 确定最大适应度对应参数
    [bestfitness, bestpop] = max(fitvalue); % 最优适应度及其对应个体位置
    % 记录每次迭代的最优，并记录对应的最优参数
    value(i, :) = BinaryToDecimal(newpop_mutation(bestpop, :), lu, num_var);
    fitness(i) = bestfitness; 
    objvalue(i) = objvalue(bestpop);
    % 更新最优
    if bestfitness >= bestfit
        bestfit = fitness(i);
        bestval = newpop_mutation(bestpop, :); % 二进制形式
        bestobjvalue = objvalue(bestpop);
    end	
    % 记录最优的变化
    bestfits(i) =  bestfit;
    % 更新种群
    pop = newpop_mutation;
end
%% 作图
figure(1)
% plot(1:1:iter_max, fitness_ave, 'r', 1:1:iter_max, y, 'b')
plot(1:1:iter_max, fitness_ave, 'r', 1:1:iter_max, fitness, 'b')
grid on
legend('每次迭代平均适应度', '每次迭代最优适应度')

figure(2)
plot(1:1:iter_max, bestfits, 'r');
grid on
legend('最优变化')

figure(3)
plot3(value(:,1), value(:,2),objvalue,'r*');                                       
hold on
[x, y1] = meshgrid(min_region:0.01:max_region);
f = @(x, y) abs( sin(pi*(x - 3))./(pi*(x - 3)) ).*abs( sin(pi*(y - 3))./(pi*(y - 3)) );
mesh(x, y1, f(x, y1));

% 输出最优结果
disp(['最优染色体为', num2str(bestval)]);
disp(['最优参数为', num2str(BinaryToDecimal(bestval, lu, num_var))]);
disp(['最优适应度为', num2str(bestfit)]);


%% 三维图

% [x, y] = meshgrid(-5:0.5:5)
% f = @(x, y) x + y + y.^(x/3) + 10.*y.*sin(5.*x.*y.*pi) + 7.*x.*y.*cos(4.*x.*pi);
% z = f(x, y)
% mesh(x, y, real(z))



%% 遗传算法仿真
% Function：
% 对一元函数求极大值(单参数)
% Created by YanBin Wang
% Date: 2024.05.16

warning off; % 关闭所有警告
clear;
clc;                    
close all

% 函数参数设置
max_region = 10;
mix_region = -10;
lu = [mix_region; max_region];
[temp, num_var] = size(lu);

f = @(x) x+10*sin(5*x)+7*cos(4*x*pi);

% 种群参数设置
popsize = 100;           	% 种群大小
chromlength = 20;          	% 二进制编码长度(个体长度)，需要根据问题求解的精度、变量的取值范围综合判定
pc = 0.8;               	% 交叉概率，只有在随机数小于pc时，才会产生交叉 一般取 60~100%
pm = 0.08;                	% 变异概率，一般 0.01~0.1
iter_max = 200;          % 迭代次数（遗传次数）

fitness_ave = -inf(iter_max, 1); % 记录每次迭代的平均适应度
bestfit = -inf;
bestfits = -inf(iter_max, 1); % 记录最优变化

% 随机生成二进制编码的初始群体，个体长度=变量的二进制编码的长度*变量个数（多参数时）           
if num_var == 1                 % 只有一个变量的情况
    pop = round(rand(popsize, chromlength)); 
else                            % 多个变量的情况
    disp('Warning：这是单参数matlab程序，请重新设置');
    exit;
end

% 迭代开始
for i = 1:1:iter_max             
    % 种群十进制形式
    pop_decimal = BinaryToDecimal(pop, lu);
    % 计算群体中每个个体的适应度
    fitvalue = calfitvalue(pop_decimal);      
    % 计算目标函数, 设定目标函数即为适应度
    objvalue = fitvalue;    
    % 选择，复制
    newpop_selection = selection(pop, fitvalue);            	
    % 交叉
    newpop_crossover = crossover(newpop_selection, pc);                
    % 变异
    newpop_mutation = mutation(newpop_crossover, pm);                 
    % 种群十进制形式
    pop_decimal = BinaryToDecimal(newpop_mutation, lu);
    % 计算群体中每个个体的适应度
    fitvalue = calfitvalue(pop_decimal);      
    % 计算目标函数, 设定目标函数即为适应度
    objvalue = fitvalue;  
    % 计算平均适应度
    fitness_ave(i) = calavefitness(fitvalue);  
    % 确定最大适应度对应参数
    [bestfitness, bestpop] = max(fitvalue); % 最优适应度及其对应个体位置
    % 记录每次迭代的最优，并记录对应的最优参数
    x(i) = BinaryToDecimal(newpop_mutation(bestpop, :), lu);
    y(i) = bestfitness; 
    % 更新最优
    if bestfitness >= bestfit
        bestfit = y(i);
        bestval = newpop_mutation(bestpop, :); % 二进制形式
    end	
    % 记录最优的变化
    bestfits(i) =  bestfit;
    % 更新种群
    pop = newpop_mutation;
end
%% 作图
figure(1)
plot(1:1:iter_max, fitness_ave, 'r', 1:1:iter_max, y, 'b')
grid on
legend('平均适应度', '最优适应度')

figure(2)
plot(1:1:iter_max, bestfits, 'r');
grid on
legend('最优变化')

figure(3)
fplot(f, lu');
hold on
plot(x, y, 'r*');                                       

%输出最优结果
disp(['最优染色体为', num2str(bestval)]);
disp(['最优参数为', num2str(BinaryToDecimal(bestval, lu))]);
disp(['最优适应度为', num2str(bestfit)]);

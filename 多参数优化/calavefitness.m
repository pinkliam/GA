%% 求平均适应度
function fitness_ave = calavefitness(fitness)
[popsize ,~] = size(fitness);
fitness_ave = sum(fitness)/popsize;
end
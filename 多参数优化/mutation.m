%% 变异
function [newpop]=mutation(pop,pm)
[popsize,chromlength]=size(pop);
newpop=ones(size(pop));
for i=1:popsize
    if(rand<pm)
        mpoint=round(rand*chromlength);              % 确定变异点
        if mpoint<=0 % 保证变异点在1--chromlength内
            mpoint=1;                   	
        end
        newpop(i,:)=pop(i,:);
        % 0/1取反变异
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


% 变异(mutation)，基因突变指父代中的每个个体的每一位都以一定概率 pm 翻转，即由“1”变“0”，或由“0”变“1”。
% 变异可以使求解过程随机地搜索到解可能存在的整个空间，可以在一定程度上求得全局最优解。
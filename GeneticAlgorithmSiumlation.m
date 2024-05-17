%% �Ŵ��㷨����
% ���ܣ� 
% һԪ������ֵ
% Created by YanBin Wang
% Date: 2024.5.16

function GeneticAlgorithmSiumlation

clear;
clc;                        % ����
close all
x_max   = 10;
x_min   = 0;
x_range = [x_min, x_max];
[num_var,temp] = size(x_range);
popsize = 100;           	% Ⱥ���С
chromlength = 20;          	% �ַ�������(���峤��)����Ҫ�����������ľ��ȡ�������ȡֵ��Χ�ۺ��ж�
pc = 0.7;               	% ������ʣ�ֻ���������С��pcʱ���Ż�������� һ��ȡ 60~100%
pm = 0.05;                	% ������ʣ�һ�� 0.1~10%
iter_num    = 200;          % �Ŵ�����

pop = initpop(popsize, chromlength, num_var);             % ������������Ʊ���ĳ�ʼȺ��
for i = 1:1:iter_num              
    % ����Ŀ�꺯��
    objvalue    = calobjvalue(pop, x_range, chromlength); 
    
    % ����Ⱥ����ÿ���������Ӧ��
    fitvalue    = calfitvalue(objvalue);             	
    
    % ѡ�񣬸���
    newpop      = selection(pop,fitvalue);            	
    
    % ����
    newpop1     = crossover(newpop,pc);                
    
    % ����
    newpop2     = mutation(newpop1,pm);                 
    
    % ����Ŀ�꺯��
    objvalue    = calobjvalue(newpop2, x_range, chromlength); 
    
    % ����Ⱥ����ÿ���������Ӧ��
    fitvalue    = calfitvalue(objvalue);            	
    
    % ƽ����Ӧ��
    fitness_ave(i) = calavefitness(fitvalue);           
    
    % ���Ⱥ������Ӧֵ���ĸ��弰����Ӧֵ
    [bestindividual,bestfit]=best(newpop2,fitvalue); 	
    
    % ���ص� y ������Ӧ��ֵ�����Ǻ���ֵ
    y(i)    = bestfit;                                              
    
    % ���Ա��������ʮ����
    x(i,:)    = decodechrom(bestindividual, x_range, chromlength);	
    
    % ������Ⱥ
    pop     = newpop2;
end

% ��ͼ
figure(1)
plot(1:iter_num, fitness_ave, 'r', 1:iter_num, y, 'b')
grid on
legend('ƽ����Ӧ��', '������Ӧ��')


figure(2)
fplot('x+10*sin(5*x)+7*cos(4*x)',x_range);
hold on
plot(x,y,'r*');                                       

%������
[z ,index]=max(y);             %�������ֵ����λ��
ymax = z;
disp(['����Ⱦɫ��Ϊ', num2str(x(index,:))])
disp(['������Ӧ��Ϊ', num2str(ymax)])

end


%% 2.1��ʼ��(����)
% initpop.m�����Ĺ�����ʵ��Ⱥ��ĳ�ʼ����popsize��ʾȺ��Ĵ�С��chromlength��ʾȾɫ��ĳ���(��ֵ���ĳ���)��
% rand�������ÿ����ԪΪ {0,1} ����Ϊpopsize������Ϊchromlength�ľ���
% round�Ծ����ÿ����Ԫ����Բ���������롣
% ���ȴ�Сȡ���ڱ����Ķ����Ʊ���ĳ���*����������
% ��ʼ��
function pop=initpop(popsize,chromlength, num_var) 
if num_var == 1                 % ֻ��һ�����������
    pop = round(rand(popsize,chromlength)); 
else                            % ������������

end

end

%% 2.2 ����Ŀ�꺯��ֵ
% 2.2.1 calobjvalue.m�����Ĺ�����ʵ��Ŀ�꺯���ļ���
% ʵ��Ŀ�꺯���ļ��㣬�� ��ֵ�� �е���ת��Ϊ ���������
% ����pop ��ʾ��Ⱥ��
% x_range ������ȡֵ��Χ
% chromlength  �����Ʊ���ĳ���
function [objvalue]=calobjvalue(pop, x_range, chromlength)
x = decodechrom(pop, x_range, chromlength);   %��popÿ��ת����ʮ������
objvalue = x+10*sin(5*x)+7*cos(4*x);   %����Ŀ�꺯��ֵ

end


% 2.2.2 �������Ʊ���ת��Ϊʮ������(2)
% decodechrom.m�����Ĺ����ǽ�Ⱦɫ��(������Ʊ���)ת��Ϊʮ���ƣ�
% ����pop ��ʾ��Ⱥ��
% x_range ������ȡֵ��Χ
% chromlength  �����Ʊ���ĳ���
% (���ڶ���������ԣ�������������������2*chromlengthλ��ʾ��ÿ������chromlengthλ
% ���һ��������1��ʼ����һ��������chromlength+1��ʼ)
function pop2=decodechrom(pop,x_range,chromlength)          %1  10

pop1 = pop;
x_min = x_range(1);
x_max = x_range(2);
pop2 = x_min + (x_max - x_min) * decodebinary(pop1)/(2^chromlength - 1);  %��popÿ��ת����ʮ����ֵ
end

% ����������ת��Ϊʮ������clc
function pop2=decodebinary(pop)
[px,py]=size(pop);                  % ��pop�к�����
for i=1:py
    pop1(:,i)=2.^(py-i).*pop(:,i);  % ���� [2^n 2^(n-1) ... 1] ����������Ȼ�����
end
pop2=sum(pop1,2);                   % ��pop1��ÿ��֮�ͣ��õ���Ӧ��ʮ������

end

%% 2.3 ����������Ӧֵ
function fitvalue=calfitvalue(objvalue)
[px,py]=size(objvalue);                   %Ŀ��ֵ�����и�
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

%% ��ƽ����Ӧ��
function fitness_ave = calavefitness(fitness)
[N ,~] = size(fitness);
fitness_ave = sum(fitness)/N;
end

%% 2.4 ѡ��
function [newpop]=selection(pop, fitvalue) 
totalfit = sum(fitvalue);                       % ����Ӧֵ֮��
fitvalue = fitvalue/totalfit;                   % �������屻ѡ��ĸ���
fitvalue = cumsum(fitvalue);                    % �� fitvalue=[1 2 3 4]����cumsum(fitvalue)=[1 3 6 10],Ҫ�ۼӣ����̶ķ������ο��Ƿ���ת�õ������� 
[row_p,col_p] = size(pop);                      % row_p*col_p
ms = sort(rand(row_p,1));                       % ��С��������
fitin = 1;
newin = 1;
while newin <= row_p                           	% ѡ��col_p���¸��壬���ظ��������������ܵķ�����̫һ��
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:) = pop(fitin,:);
        newin = newin + 1;
    else
        fitin = fitin + 1;
    end
end

end


%% 2.5 ����
function [newpop]=crossover(pop,pc)                         % pc=0.6
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:2:px-1                                              % ����Ϊ2���ǽ����ڵ�����������н���
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

%% 2.6 ����
function [newpop]=mutation(pop,pm)
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:px
    if(rand<pm)
        mpoint=round(rand*py);              % �����ı������1-10֮��
        if mpoint<=0
            mpoint=1;                   	% ����λ��
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


%% 2.7 ���Ⱥ����������Ӧֵ�������
% �Ŵ��㷨�ӳ���
% ����� t ��Ⱥ������Ӧֵ����ֵ
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



function fitvalue=calfitvalue(pop_decimal)
% 计算个体的适应值
f = @(x) x + 10*sin(5*x) + 7*cos(4*x*pi);

fitvalue = f(pop_decimal);

end
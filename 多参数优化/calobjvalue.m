%% 计算目标函数
function objvalue=calobjvalue(pop_decimal)

[popsize, ~] = size(pop_decimal);
f = @(x, y) abs( sin(pi*(x - 3))./(pi*(x - 3)) ).*abs( sin(pi*(y - 3))./(pi*(y - 3)) );

objvalue= f(pop_decimal(:, 1), pop_decimal(:, 2));

end
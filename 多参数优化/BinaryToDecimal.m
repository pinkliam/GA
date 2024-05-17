function pop_decimal=BinaryToDecimal(pop,lu, num_var)  
% 将种群二进制转换成十进制
[popsize, chromlength] = size(pop);
temp_pop1 = zeros(popsize, chromlength);
min_region = lu(1, :);
max_region = lu(2, :);

% 将pop每行的每个变量均转换成十进制值
for i = 1:1:chromlength
    k = mod(i, chromlength/num_var); % 取余数
    if k == 0
        k = next_k;  % 没有余数时（对应此个变量的最后一位二进制）
    end
    next_k = k + 1;
    temp_pop1(:,i) = 2^(k-1).*pop(:, i);       
end
for j = 1:1:num_var
    temp_pop2(:, j) = sum(temp_pop1(:, ((chromlength*(j-1))/num_var+1):((chromlength*j)/num_var)), 2);
end
pop_decimal = temp_pop2.*repmat((max_region - min_region)./(2^(chromlength/num_var) - 1), popsize,1) + repmat(lu(1, :), popsize,1);

end
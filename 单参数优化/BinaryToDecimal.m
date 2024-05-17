function pop_decimal=BinaryToDecimal(pop,lu)  
% 将种群二进制转换成十进制
[popsize, chromlength] = size(pop);
temp = zeros(popsize, 1);
min_region = lu(1);
max_region = lu(2);

% 将pop每行转换成十进制值
for i = 1:1:chromlength
    temp = temp + 2^(i-1).*pop(:, i);
end
pop_decimal = temp.*(max_region - min_region)./(2^chromlength - 1) + lu(1); % 防止是负数，要加上最小的负数边界

end
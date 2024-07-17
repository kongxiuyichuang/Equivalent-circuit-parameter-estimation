function fit = seri_evalFitness(pop,Obj_function)
% pop: 当前粒子群的位置矩阵 (popsize x nvars)

% 粒子群大小和变量维度
[popsize, nvars] = size(pop);

% 初始化适应度函数值向量
fit = zeros(popsize, 1);

% 计算适应度函数值
for i = 1:popsize
    fun=Obj_function;
    % 将每个变量作为参数，调用适应度函数计算适应度函数值
    x = pop(i,:);
    fit(i) =fun(x);
end

end


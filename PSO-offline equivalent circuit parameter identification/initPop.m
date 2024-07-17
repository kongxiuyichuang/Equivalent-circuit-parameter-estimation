function pop = initPop(popsize, nvars, lb, ub)
% popsize: 粒子群大小
% nvars: 变量维度
% lb: 变量下限向量 (1 x nvars)
% ub: 变量上限向量 (1 x nvars)

% 随机数生成器
rng('shuffle');

% 初始化粒子群的位置
pop = rand(popsize, nvars) .* repmat(ub-lb, popsize, 1) + repmat(lb, popsize, 1);

end

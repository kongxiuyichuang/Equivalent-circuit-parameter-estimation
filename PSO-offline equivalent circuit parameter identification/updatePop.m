function [pop, vel] = updatePop(pop, vel, pbest, gbest, w, c1, c2, lb, ub, VelMin, VelMax)
% pop: 当前粒子群的位置矩阵 (popsize x nvars)
% vel: 当前粒子群的速度矩阵 (popsize x nvars)
% pbest: 当前粒子群的个体最优解矩阵 (popsize x nvars)
% gbest: 当前粒子群的全局最优解向量 (1 x nvars)
% w: 惯性权重因子
% c1: 自身经验因子
% c2: 社会经验因子
% lb: 变量下限向量 (1 x nvars)
% ub: 变量上限向量 (1 x nvars)

% 粒子群大小和变量维度
[popsize, nvars] = size(pop);

% 随机数生成器
rng('shuffle');

% 更新速度和位置
for i = 1:popsize
    % 计算惯性项，自身经验项和社会经验项
    r1 = rand(1, nvars);
    r2 = rand(1, nvars);
    vel(i,:) = w * vel(i,:) ...
        + c1 * r1 .* (pbest(i,:) - pop(i,:)) ...
        + c2 * r2 .* (gbest - pop(i,:));

    % 控制速度和位置越界
    % vel(i,:) = max(min(vel(i,:), ub - lb), lb - ub);
    % Apply Velocity Limits
    vel(i,:) = max( vel(i,:),VelMin);
    vel(i,:) = min( vel(i,:),VelMax);

    % Update Position
    pop(i,:) = pop(i,:) +  vel(i,:);

    % Velocity Mirror Effect
    IsOutside=(pop(i,:)<lb | pop(i,:)>ub);
    vel(IsOutside,:)=- vel(IsOutside,:);

    % Apply Position Limits
    pop(i,:) = max(pop(i,:),lb);
    pop(i,:) = min(pop(i,:),ub);

    %pop(i,:) = max(min(pop(i,:) + vel(i,:), ub), lb);
end

end

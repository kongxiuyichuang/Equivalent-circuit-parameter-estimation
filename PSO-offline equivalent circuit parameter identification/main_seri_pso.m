function [history,gbest,gbestcost]=main_seri_pso(popsize,maxiter,lb,ub,Obj_function)
%seri-pso main function
% -----Input-------
% popsize:swarm size
% maxiter:max number of iteration
% lb:lower bound of decision variable
% ub:upper bound of decision variable
% Obj_function:cost function
% -----Output--------
% history:optimization record
% gbest:best solution
% gbestcost: best cost

tic
%% Intialize
nvars=length(lb);
pop = initPop(popsize, nvars, lb, ub);
% 计算适应度函数值
fit = seri_evalFitness(pop,Obj_function);
% 初始化个体最优解向量和全局最优解向量
pbest = pop;
pbestcost=fit;
[gbestcost,gbest_idx]=min(fit);
gbest=pop(gbest_idx,:);
bestSol_id_iter_number=0;
bestSol_id_particle_number=gbest_idx;
%粒子的速度进行初始化
% vel = rand(popsize, nvars) .* repmat((ub-lb)/2, popsize, 1);
vel =zeros(popsize,nvars);

% 初始化速度和适应度函数值结构体
history.vel = cell(maxiter, 1);
history.pos = cell(maxiter, 1);
history.fitness = cell(maxiter, 1);
history.bestcost=zeros(maxiter,1);
history.bestposition=zeros(maxiter,nvars);
history.bestSol_id_iter_number=cell(maxiter, 1);
history.bestSol_id_particle_number=cell(maxiter, 1);

%Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wdamp=1;        % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient
% Velocity Limits
VelMax=0.1*(ub-lb);
VelMin=-VelMax;

%% PSO Main Loop
for i = 1:maxiter
    % Updata Position and Velocity
    [pop, vel] = updatePop(pop, vel, pbest, gbest, w, c1, c2,lb, ub, VelMin, VelMax);

    % Evaluation object function
    fit = seri_evalFitness(pop,Obj_function);

    % Updata Personal best
    for j = 1:popsize      
        if   fit(j) < pbestcost(j)
            pbest(j,:) = pop(j,:);
        end
    end

   % Updata Global Best
    [minfit, idx] = min(fit);
    if minfit < gbestcost
        gbest = pop(idx,:);
        gbestcost=minfit;
        bestSol_id_iter_number=i;
        bestSol_id_particle_number=idx;
    end
    w=w*wdamp;

    % Record speed\position\fitness function values\bestposition\bestcost
    history.vel{i} = vel;
    history.pos{i} = pop;
    history.fitness{i} = fit;

    history.bestposition(i,:)=gbest;
    history.bestcost(i)=gbestcost;
    history.bestSol_id_iter_number{i}=bestSol_id_iter_number;
    history.bestSol_id_particle_number{i}=bestSol_id_particle_number;

    % Display current progress
    fprintf('Iteration %d: Best fitness = %f\n', i, history.bestcost(i));

end
%% plot 
% figure;
% semilogy( history.bestcost,'-o','LineWidth',2);
% xlabel('Iteration');ylabel('Best Cost');grid on;
% save('seri_pso_recoder.mat');
toc
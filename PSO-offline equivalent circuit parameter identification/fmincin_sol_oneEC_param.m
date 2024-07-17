%% load SOC-OCV data
clear
close all
clc
%% input hppc_data
load("SOC_OCV.mat")
hppc_data=readmatrix('hppc_p&n_raw_data.txt');

[~,ia,~]=unique(hppc_data(:,1),'stable');
soc=hppc_data(ia,1);
data=struct();
for i=1:length(ia)
    if i==length(ia)
        data(i).soc=soc(i);
        data(i).time=hppc_data(ia(i):length(hppc_data),2);
        data(i).current=hppc_data(ia(i):length(hppc_data),3);
        data(i).voltage=hppc_data(ia(i):length(hppc_data),4);
    else
        data(i).soc=soc(i);
        data(i).time=hppc_data(ia(i):ia(i+1)-1,2);
        data(i).current=hppc_data(ia(i):ia(i+1)-1,3);
        data(i).voltage=hppc_data(ia(i):ia(i+1)-1,4);
    end
end
delta_t=hppc_data(2,2)-hppc_data(1,2);

%% solution EC parameters
sol=struct();
fmincon_result=struct();
history=cell(length(data),1);
gbest=cell(length(data),1);
gbestcost=cell(length(data),1);
order=1;
filename=["fmincon_result","pso_result"];

for i=1:length(data)
    SOC_0=data(i).soc;
    current=data(i).current;
    t=data(i).time;
    Ut_exp=data(i).voltage;
    fun=@(val_sol) object_fun(SOC_OCV,SOC_0,current,val_sol,t,delta_t,Ut_exp,order);

    % fmincon求解
    ub=[0.003,1e6,0.01];
    lb=[0.0001,1e3,1e-4];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    x_0=[0.001,1e4,0.001];
    [sol(i).x,sol(i).fval] = fmincon(fun,x_0,A,b,Aeq,beq,lb,ub);

    % pso求解
    VarMax=[0.003,1e6,0.01];
    VarMin=[0.0001,1e3,1e-4];
    nVar=length(VarMin);
    MaxIt=100;
    nPop=30;
    
    [history{i},gbest{i},gbestcost{i}]=main_seri_pso(nPop,MaxIt,VarMin,VarMax,fun);

    figure()
    plot(t,Ut_exp,'--ob',DisplayName='exp');hold on;
    %plot pso result
    plot(t,frist_order_EC_model(SOC_OCV,SOC_0,current,...
        gbest{i}(1),gbest{i}(2),gbest{i}(3),t,delta_t),'-r',LineWidth=2,DisplayName='pso-sim');hold on;
    %plot fmincon result
    plot(t,frist_order_EC_model(SOC_OCV,SOC_0,current,...
        sol(i).x(1),sol(i).x(2),sol(i).x(3),t,delta_t),'-g',LineWidth=2,DisplayName='fimcon-sim');

    grid on;legend(Location="southeast");xlabel('time(s)');ylabel('voltage(V)')
    label=sprintf('soc=%s%%',num2str(SOC_0*100));title(label);
end

% fmincon sol param initialization
fmincon_result_R_0=zeros(length(sol),1);
fmincon_result_C_p=zeros(length(sol),1);
fmincon_result_R_p=zeros(length(sol),1);
fmincon_result_SOC=zeros(length(sol),1);
fmincon_result_tau=zeros(length(sol),1);
% pso sol param initialization
R_0=zeros(length(gbest),1);
C_p=zeros(length(gbest),1);
R_p=zeros(length(gbest),1);
SOC=zeros(length(gbest),1);
tau=zeros(length(gbest),1);
for i=1:length(sol)
    fmincon_result_R_0(i)=sol(i).x(1);
    fmincon_result_C_p(i)=sol(i).x(2);
    fmincon_result_R_p(i)=sol(i).x(3);
    fmincon_result_SOC(i)=data(i).soc;
    fmincon_result_tau(i)=fmincon_result_C_p(i)*fmincon_result_R_p(i);

    R_0(i)=gbest{i}(1);
    C_p(i)=gbest{i}(2);
    R_p(i)=gbest{i}(3);
    SOC(i)=data(i).soc;
    tau(i)=C_p(i)*R_p(i);
end
%% output EC_parameter
fmincon_EC_para=[fmincon_result_SOC,fmincon_result_R_0,fmincon_result_C_p,fmincon_result_R_p];
EC_para=[SOC,R_0,C_p,R_p];
%save('EC_para.mat',"EC_para");

%% plot parameter result
figure();
plot(SOC,R_0,'-o',DisplayName='pso-R_0'); hold on 
plot(fmincon_result_SOC,fmincon_result_R_0,'-o',DisplayName='fmincon-R_0'); hold on 
grid on;legend;xlabel('SOC');ylabel('R_0(\Omega)') ;

figure();
plot(SOC,C_p,'-o',DisplayName='pso-C_p');hold on
plot(fmincon_result_SOC,fmincon_result_C_p,'-o',DisplayName='fmincon-C_p');
grid on;legend;xlabel('SOC');ylabel('C_p(F)') ;

figure();
plot(SOC,R_p,'-o',DisplayName='pso-R_p');hold on;
plot(fmincon_result_SOC,fmincon_result_R_p,'-o',DisplayName='fminconR_p')
grid on;legend;xlabel('SOC');ylabel('R_p(\Omega)') ;

figure();
plot(SOC,tau,'-o',DisplayName='pso-tau');hold on
plot(fmincon_result_SOC,fmincon_result_tau,'-o',DisplayName='fmincon-tau')
grid on;legend;xlabel('SOC');ylabel('tau(s)') ;

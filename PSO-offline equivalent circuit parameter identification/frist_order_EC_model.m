function U_t=frist_order_EC_model(SOC_OCV,SOC_0,current,R_0,C_p,R_p,t,delta_t)
%一阶等效电路模型,假设脉冲期间电池SOC变化忽略不变
% input:
% SOC_OCV:SOC OCV矩阵，size:N*2
% SOC_0:起始SOC
% current:负载电流
% R_0：欧姆内阻
% C_p：极化电容
% R_p：极化电阻
% delta_t：时间离散间隔
% t：时间序列从0开始，size:N*1
% output:
% U_t:终端电压,size:N*1

tau=R_p*C_p;
OCV=interp1(SOC_OCV(:,1),SOC_OCV(:,2),SOC_0,"linear","extrap");

U_p=zeros(length(t),1);
U_t=zeros(length(t),1);
U_p(1)=0;
U_t(1)=OCV;
for k=2:length(t)
    U_p(k)=U_p(k-1)*exp(-delta_t/tau)+R_p*(current(k-1))*(1-exp(-delta_t/tau));
    U_t(k)=OCV-U_p(k)-current(k)*R_0;

end

end
function U_t=second_order_EC_model(SOC_OCV,SOC_0,current,R_0,C_p1,R_p1,C_p2,R_p2,t,delta_t)
%二阶等效电路模型
% t,时间序列从0开始
tau1=R_p1*C_p1;
tau2=R_p2*C_p2;
OCV=interp1(SOC_OCV(:,1),SOC_OCV(:,2),SOC_0,"linear","extrap");
%delta_t=1;
U_p1=zeros(length(t),1);
U_p2=zeros(length(t),1);
U_t=zeros(length(t),1);
U_p1(1)=0;
U_p2(1)=0;
U_t(1)=OCV;
for k=2:length(t)
U_p1(k)=U_p1(k-1)*exp(-delta_t/tau1)+R_p1*(current(k-1))*(1-exp(-delta_t/tau1));
U_p2(k)=U_p2(k-1)*exp(-delta_t/tau2)+R_p2*(current(k-1))*(1-exp(-delta_t/tau2));

U_t(k)=OCV-U_p1(k)-U_p2(k)-current(k)*R_0;
end
end
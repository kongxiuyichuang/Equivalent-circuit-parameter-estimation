function RMSE=object_fun(SOC_OCV,SOC_0,current,val_sol,t,delta_t,Ut_exp,order)
% 等效电路参数辨识目标函数
%input：
% SOC_OCV:SOC OCV矩阵，size:N*2
% SOC_0:起始SOC
% current:负载电流
% val_sol：1阶时,size:1*3,[R_0,C_p,R_p];2阶时,size:1*5,[R_0,C_p1,R_p1,C_p2,R_p2]
% t：时间序列从0开始，size:N*1
% delta_t：时间离散间隔
% Ut_exp:实测端电压序列，size:N*1
% order: 阶数,1 & 2
% output:
% RMSE:模型计算值与实验值电压均方差

if order==1
    R_0=val_sol(1);
    C_p=val_sol(2);
    R_p=val_sol(3);
    Ut_sim=frist_order_EC_model(SOC_OCV,SOC_0,current,R_0,C_p,R_p,t,delta_t);
elseif order==2
    R_0=val_sol(1);
    C_p1=val_sol(2);
    R_p1=val_sol(3);
    C_p2=val_sol(4);
    R_p2=val_sol(5);
    Ut_sim=second_order_EC_model(SOC_OCV,SOC_0,current,R_0,C_p1,R_p1,C_p2,R_p2,t,delta_t);
else
    error('error: order is not 1 or 2')
end
diff=Ut_exp-Ut_sim;
RMSE=sqrt(mean(diff.*diff));
end
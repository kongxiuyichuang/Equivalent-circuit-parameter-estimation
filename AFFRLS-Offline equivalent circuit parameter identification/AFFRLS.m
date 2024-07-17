%y 是输出数据Ut，x 是输入数据I，lambda 是遗忘因子，theta_init 是初始参数向量。这个函数返回估计的参数向量 theta
%theta=[(1-alph_1)*U_OCV,alph_1,alph_2,alph_3]'
%Phil=[1,U_t_{k-1},I_L_{k},I_L_{k-1}]
%Ut=Phil*theta;

clear
close all

%% load hppc data
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
        data(i).OCV=hppc_data(ia(i),4);
    else
        data(i).soc=soc(i);
        data(i).time=hppc_data(ia(i):ia(i+1)-1,2);
        data(i).current=hppc_data(ia(i):ia(i+1)-1,3);
        data(i).voltage=hppc_data(ia(i):ia(i+1)-1,4);
        data(i).OCV=hppc_data(ia(i),4);
    end
end

%% RLS Process
RLS_result=struct();
OCV=zeros(length(data),1);
OCV_exp=zeros(length(data),1);
voltage_erro=zeros(length(data),1);

figure()
for k=1:length(data)

    y=data(k).voltage;
    I=data(k).current;
    time=data(k).time;
    %soc=data(k).soc;
    OCV_exp(k)=data(k).OCV;
    %% RLS parameter init
    theta_init=zeros(4,1);
    Phil=zeros(1,4);
    lambda=0.97;
    K=zeros(4,1);
    P = 10^6*eye(length(theta_init));
    N = length(y);

    %% output parameter init
    U_oc=zeros(N,1);
    Rs=zeros(N,1);
    tau_1=zeros(N,1);
    Rp=zeros(N,1);
    Cp=zeros(N,1);
    Voltage_erro=zeros(N,1);
    y_pre=zeros(N,1);
    e=zeros(N,1);
    theta=theta_init;
    for t = 2:N
        Phil=[1,y(t-1),I(t),I(t-1)];
        % 计算预测误差
        y_pre(t)=Phil*theta;
        e(t) = y(t) - y_pre(t) ;
        %自适应遗忘因子
        if abs(I(t)-I(t-1))>2*min(abs(I(t)),abs(I(t-1)))
            lambda=1-e(t-1).^2/(1+ Phil* P * Phil');
        else
            lambda=1;
        end
        % 计算增益矩阵
        K = P * Phil' / (lambda + Phil* P * Phil');
        % 更新参数向量
        theta = theta + K * e(t);
        % 更新协方差矩阵
        P = (P - K * Phil * P) / lambda;
        %辨识参数解析
        U_oc(t)=theta(1)/(1-theta(2));
        Voltage_erro(t)=e(t);
        Rp(t)=(2*(theta(4)+theta(2)*theta(3)))/(theta(2)^2-1);
        Rs(t)=(theta(4)-theta(3))/(theta(2)+1);
        Cp(t)=-((time(2)-time(1))*theta(2)^2+2*(time(2)-time(1))*theta(2)+(time(2)-time(1)))/(4*theta(4)+4*theta(2)*theta(3));
        tau_1(t)=Cp(t)*Rp(t);
    end
    RLS_result(k).OCV=mean(U_oc(3:end));
    RLS_result(k).Rs=Rs(end);
    RLS_result(k).Rp=Rp(end);
    RLS_result(k).Cp=Cp(end);
    RLS_result(k).tau=tau_1(end);
    % RLS_result(k).Rs_mean=mean(Rs);
    % RLS_result(k).Rp_mean=mean(Rp);
    % RLS_result(k).Cp_mean=mean(Cp);
    % RLS_result(k).tau_mean=mean(tau_1);
    RLS_result(k).voltage_RMSE=sqrt(mean(Voltage_erro(3:end).^2));
    RLS_result(k).OCV_exp=OCV_exp(k);
    OCV(k)=RLS_result(k).OCV;

    %% plot
    subplot(2,5,k)
    plot(time(3:end),y(3:end),'-ob',DisplayName='Ut-exp');hold on;
    plot(time(3:end),y_pre(3:end),'-r',LineWidth=2,DisplayName='Ut-sim');
    plot(time(3:end),U_oc(3:end),'-g',LineWidth=2,DisplayName='OCV-sim');
    %legend("Location","southeast")
    label=sprintf('soc=%s%%',num2str(soc(k)*100));title(label);
    xlabel('time(s)');ylabel('Voltage(V)')
end
legend("Location","southeast")

figure()
yyaxis left
plot(soc,OCV,'o-b',DisplayName='RLS-result');hold on;
plot(soc,OCV_exp,'-*r',DisplayName='exp-result');
xlabel('SOC');ylabel('OCV');
yyaxis right
plot(soc,(OCV_exp-OCV)*1000,'-^',DisplayName='OCV-erro');
ylabel('OCV-erro(mV)');legend("Location","northwest");grid on;
save('AFFRLS_result.mat','RLS_result')

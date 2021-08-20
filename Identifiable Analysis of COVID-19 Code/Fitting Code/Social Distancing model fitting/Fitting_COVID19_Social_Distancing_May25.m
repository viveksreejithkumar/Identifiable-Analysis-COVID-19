clear all
close all
clc

global  epiData COVID_Deaths epiForward epiMeasure  initial_cond N

epiData = [193;300;516;561;693;725;917;899;945;1100;1300;1300;1200;764;...
          1200;1000;1100;1000;1100;1100;814;925;593;903;1200;1100;743;...
          713;749;860;706;1200;739;694;624;593;341;514;978;709;632;770;...
          575;583;626;343;768;683;384;608;802;794;790;482;1142;661;615;...
          477;1209;739;753;563];

COVID_Deaths = [5;13;20;13;19;12;19;22;31;31;35;30;36;41;44;43;41;37;42;...
41;35;45;49;43;43;56;39;49;49;47;43;52;39;41;50;42;51;44;51;52;30;40;59;...
38;49;41;34;32;38;41;44;38;46;33;33;33;22;33;34;29;24;35];

tdata = 0:1:61;

dt = 0.1; tf = 61;

epiForward = (0:dt:tf);
epiMeasure = 1:1/dt:length(epiForward);
N = 210000;
initial_cond = [N-193 0 193 0 0];
initial_cond_ns = [N-193 193 0 0];


%k = [beta eta delta gamma nu rho betaE deltaE k];
% k = [0.21014 0.6 0.01 0.3 0.009 1/60 0.2];
% k = [1.20724887779556,0.199215122470033,0.000344641632357087,...
%     0.288460308248421,0.00627214924740963,0.0830732791660550,0.603123799509962];
% k = [0.668480531160839*1.66941701350445,0.130006113573208,0.0189844324736501,...
%     0.111114633491532,0.00829554667499385,0.0247452661226331,...
%     0.668480531160839,0.0189844324736501,1/14]

params = [0.0423722646592824,0.156181646449812,0.0229258085719113,...
    0.178586841700095,0.0110252925203845,0.0172801497190386,...
    0.787270389284732,0.0437398345879143,0.0948158165896732];
%results
%params =
%[0.00398555771219389,0.398515014117417,0.148671772481456,0.199895044751912,0.00798592714501395,0.0529338104456878,1.59638257459127,0.161583566704927,0.392442624261959];
lb = zeros(size(params));
  
lb = [0 0 0 0.0 0 0 0 0 0.0] 
ub = [1e+4 1 1 0.2 1 1e+3 1e+4 1 0.5] 
%[params,fval] =  fminsearchbnd(@err_in_data,params,lb,ub,optimset('Display','iter',...
%   'MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-8,'TolX',1e-8))


[~, y_r] = ode15s(@(t,y)Model_SIR_ODE(y,params),epiForward,initial_cond);
[~, y_no_soc_dis] = ode15s(@(t,y)Model_SIR_NO_SOCIAL_DISTANCING(y,params),epiForward,initial_cond_ns);

figure(1)
plot(epiForward,params(9)*y_r(:,3),'LineWidth',2)
hold on 
plot(tdata, epiData, 'r.', 'MarkerSize',20)
title({'COVID-19 Incidences with Social Distancing', 'March 23-May 25'}, 'Fontsize',14)
xlabel('Days since March 23','Fontsize',14)
ylabel('Incidences','Fontsize',14)
set(gca,'linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',14)
xline(7,'--k',{'March 30th, 2020'})
xline(14,'--k',{'April 6th, 2020'})
xline(21,'--k',{'April 13th, 2020'})
xline(28,'--k',{'April 20th, 2020'})
xline(35,'--k',{'April 27th, 2020'})
xline(42,'--k',{'May 4th, 2020'})
xline(49,'--k',{'May 11th, 2020'})
xline(56,'--k',{'May 18th, 2020'})
xline(63,'--k',{'May 25th, 2020'})
xticks([ 7 14 21 28 35 42 49 56 63])
xticklabels({'3/30',' 4/6','4/13', '4/20','4/27','5/4','5/11','5/18','5/25'})

figure(2)
plot(epiForward,params(5)*y_r(:,4),'LineWidth',2.5)
hold on 
plot(tdata, COVID_Deaths, 'r.', 'MarkerSize',20)
title({'COVID-19 Deaths with Social Distancing', 'March 23-May 25'}, 'Fontsize',14)
xlabel('Days since March 23','Fontsize',14)
ylabel('Deaths','Fontsize',14)
set(gca,'linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',14)
xline(7,'--k',{'March 30th, 2020'})
xline(14,'--k',{'April 6th, 2020'})
xline(21,'--k',{'April 13th, 2020'})
xline(28,'--k',{'April 20th, 2020'})
xline(35,'--k',{'April 27th, 2020'})
xline(42,'--k',{'May 4th, 2020'})
xline(49,'--k',{'May 11th, 2020'})
xline(56,'--k',{'May 18th, 2020'})
xline(63,'--k',{'May 25th, 2020'})
xticks([ 7 14 21 28 35 42 49 56 63])
xticklabels({'3/30',' 4/6','4/13', '4/20','4/27','5/4','5/11','5/18','5/25'})

figure(3)
plot(epiForward,params(9)*y_no_soc_dis(:,2),'LineWidth',2)
hold on 
plot(tdata, epiData, 'r.', 'MarkerSize',20)
title({'COVID-19 Incidences without Social Distancing', 'March 23-May 25'}, 'Fontsize',14)
xlabel('Days since March 23','Fontsize',14)
ylabel('Incidences','Fontsize',14)
set(gca,'linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',14)
xline(7,'--k',{'March 30th, 2020'})
xline(14,'--k',{'April 6th, 2020'})
xline(21,'--k',{'April 13th, 2020'})
xline(28,'--k',{'April 20th, 2020'})
xline(35,'--k',{'April 27th, 2020'})
xline(42,'--k',{'May 4th, 2020'})
xline(49,'--k',{'May 11th, 2020'})
xline(56,'--k',{'May 18th, 2020'})
xline(63,'--k',{'May 25th, 2020'})
xticks([ 7 14 21 28 35 42 49 56 63])
xticklabels({'3/30',' 4/6','4/13', '4/20','4/27','5/4','5/11','5/18','5/25'})

figure(4)
plot(epiForward,params(5)*y_no_soc_dis(:,3),'LineWidth',2.5)
hold on 
plot(tdata, COVID_Deaths, 'r.', 'MarkerSize',20)
title({'COVID-19 Deaths without Social Distancing', 'March 23-May 25'}, 'Fontsize',14)
xlabel('Days since March 23','Fontsize',14)
ylabel('Deaths','Fontsize',14)
set(gca,'linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',14)
xline(7,'--k',{'March 30th, 2020'})
xline(14,'--k',{'April 6th, 2020'})
xline(21,'--k',{'April 13th, 2020'})
xline(28,'--k',{'April 20th, 2020'})
xline(35,'--k',{'April 27th, 2020'})
xline(42,'--k',{'May 4th, 2020'})
xline(49,'--k',{'May 11th, 2020'})
xline(56,'--k',{'May 18th, 2020'})
xline(63,'--k',{'May 25th, 2020'})
xticks([ 7 14 21 28 35 42 49 56 63])
xticklabels({'3/30',' 4/6','4/13', '4/20','4/27','5/4','5/11','5/18','5/25'})

function error_in_data = err_in_data(k) 

global  epiData COVID_Deaths epiMeasure epiForward initial_cond N

 [~,y] = ode15s(@(t,y)Model_SIR_ODE(y,k),epiForward,initial_cond);
 

 
 Incidences = k(9)*y(epiMeasure(:),3);
 
 Deaths = k(5)*y(epiMeasure(:),4);                
 error_in_data = (Incidences - epiData)'*(Incidences - epiData)+...
                 ((Deaths - COVID_Deaths)'*(Deaths - COVID_Deaths))*100;           
  
end

function dy = Model_SIR_ODE(y,k)
global N


dy = zeros(5,1);

beta = k(1);
eta = k(2);
delta = k(3);
gamma = k(4);
nu = k(5);
rho = k(6);
betaE = k(7);
deltaE = k(8);
k = k(9);

S = y(1);
Sd = y(2);
E = y(3);
I = y(4);
R = y(5);


dy(1) =  - betaE*S.*E./N - beta*S.*I./N - eta*S + rho*Sd;
dy(2) = eta*S  - betaE*deltaE*Sd.*E./N - beta*delta*Sd.*I./N - rho*Sd;
dy(3) = beta*(S + delta*Sd).*I./N + betaE*(S + deltaE*Sd).*E./N - k*E;
dy(4) = k*E - (gamma + nu)*I;
dy(5) = gamma*I; 
end

function dy = Model_SIR_NO_SOCIAL_DISTANCING(y,k)
global N


dy = zeros(4,1);

beta = k(1);
eta = k(2);
delta = k(3);
gamma = k(4);
nu = k(5);
rho = k(6);
betaE = k(7);
deltaE = k(8);
k = k(9);

S = y(1);
E = y(2);
I = y(3);
R = y(4);


dy(1) =  - betaE*S.*E./N - beta*S.*I./N;
dy(2) = betaE*S.*E./N + beta*S.*I./N - k*E;
dy(3) = k*E - (gamma + nu)*I;
dy(4) = gamma*I;

end
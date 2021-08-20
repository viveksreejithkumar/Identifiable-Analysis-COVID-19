clear all
close all
clc

global  epiData epiForward 

epiData = [193
300
516
561
693
725
917
899
945
1100
1300
1300
1200
764
1200
1000
1100
1000
1100
1100
814
925
593
903
1200
1100
743
713
749
860
706
1200
739
694
624
593
341
514
978
709
631
769
575
582
626
344
766
681
384
606
800
794
787
482
1100
655
613
475
1200
727
733
545
981
503
477
614
1200
753
870
702
602
1208
1265
1020
1246
1391
1026
1147
1232
1669
1640
2322
2124
1767
2693
2416
2696
3316
4007
4628
2756
3822
5481
5013
8849
9556
8370
5378
6035
6553
9459
9373
11370
9896
6257
7301
9985
8924
11293
10316
15184
12232
9216
10015
13789
11224
10207
12418
10430
9282
9676
10129
12353
12075
9269
8828
9164
9352
9873
8927
9547
7080
4822
5382
5438
7603
7576
8403
6178
4160
5802
8123
6252
6151
6312
3829
2704
3875
4150
4546
4640
4313
3007
2051
2736
3266
3174
3581
3112
2610
1682
7420
2499
3589
3216
3699
2473
1870
1802
2324
2535
3652
3160
2376
1779
2650
2300
3267
3187
3556
2499
1656
2391
2590
2593
2809
2766
1841
808
3184
2052
2798
2664
2772
1857
1379
2261
2501
3258
2866
1779
3629
1563
2679
2835
3287
3380
3977
2459
1740
3580
2137
5427
3692
4388
2341
3322
4253
4038
4152
5406
2498
6486
2668
4561
4356
6128
5162
4359
6620
3832
4280
5685
5488
6753
4428
9896
4519
7312
7737
8955
8844
8188
6483
6093
8091
8104
10144
6621
6133
7189
6466
8617
9892
10638
9984
10190
8207
7529
7824
9376
11240
11356
10415
8854
8308
9419
11312
12880
12941
11428
8187
10717
10204
11113
12858
10653
5954
7060
8026
11879
13449
16592
20976
9737
10372
11215];


tdata = 0:1:286;

dt = 0.1; tf = 286;

epiForward = (0:dt:tf);

epiForward_vaccine = tf:dt:tf+210;

params = [0.0288252666791333,0.199999988628206,0.00206965866111908,0.524044990857231,0.499999992324709];

%create the spline from incidence data 
pp = spline(tdata,epiData);
%find the derivative
p_der=fnder(pp,1);

%try smoothing the data

smooth_epiData = smoothdata(epiData);

%create the spline from smoothincidene data 
pp_smooth = spline(tdata,smooth_epiData);
%find the derivative
p_der_smooth=fnder(pp_smooth,1);
p_secder_smooth=fnder(pp_smooth,2);
time_beta = (0:dt:286);
beta_t =  (ppval(p_der_smooth,time_beta) + params(5)*ppval(pp_smooth,time_beta))/params(5);
%beta_t =  (ppval(p_secder_smooth,time_beta) + (1+params(5))*ppval(p_der_smooth,time_beta) + params(5)*ppval(pp_smooth,time_beta))/params(5);
final_beta_t = beta_t(end);

N = 21000000;%650000*4.5;
initial_cond = [N-23 23 1300 32000];

[~, y_r] = ode15s(@(t,y)Model_SIR_ODE(t,y,params,epiForward,beta_t),epiForward,initial_cond);

new_initial = [y_r(end,1) 0 y_r(end,2)  y_r(end,3)  y_r(end,4)];
constant = y_r(end,1)*(y_r(end,3)/N + (params(4)/params(5))*y_r(end,2)/N) 

[~, y_v] = ode15s(@(t,y)Model_SIR_ODE_vaccination(y,params,final_beta_t),epiForward_vaccine,new_initial);

dimc = [0.3 0.3 0.3];

figure(1)
plot(epiForward,params(5)*y_r(:,2),'LineWidth',2)
hold on 
bar(tdata, epiData,'FaceColor',dimc)
%plot(tdata, epiData, 'r.', 'MarkerSize',20)
title({'COVID-19 Incidences (March 23-January 3)'},'fontweight','normal','fontsize',20)
xlabel('Dates starting from March 23','fontweight','normal','fontsize',15)
ylabel('Number of Daily COVID-19 Incidences','fontweight','normal','fontsize',15)
xline(7,'--k',{'April 2020'})
xline(37,'--k',{'May 2020'})
xline(67,'--k',{'June 2020'})
xline(97,'--k',{'July 2020'})
xline(127,'--k',{'August 2020'})
xline(157,'--k',{'Sept 2020'})
xline(187,'--k',{'Oct 2020'})
xline(217,'--k',{'Nov 2020'})
xline(247,'--k',{'Dec 2020'})
xline(277,'--k',{'Jan 2021'})
xticks([ 7 37 67 97 127 157 187 217 247 277])
xticklabels({'4/2020',' 5/2020','6/2020', '7/2020','8/2020','9/2020','10/2020','11/2020','12/2020','1/2021'})
yticks([ 3000 6000 9000 12000 15000 18000 21000])
yticklabels({'3000', '6000', '9000', '12000', '15000','18000','21000'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')




figure(2)
plot(tdata, epiData,'o',epiForward,ppval(pp,epiForward))
hold on 
title({'B-Spline Interpolation of Florida Department of Health Data (Without Smoothing)'}, 'Fontsize',20)
% hold on 
% plot(epiForward,ppval(p_der,epiForward),'-k')
xlabel('Dates starting from March 23','fontweight','normal','fontsize',15)
ylabel('Number of Daily COVID-19 Incidences','fontweight','normal','fontsize',15)
xline(7,'--k',{'April 2020'})
xline(37,'--k',{'May 2020'})
xline(67,'--k',{'June 2020'})
xline(97,'--k',{'July 2020'})
xline(127,'--k',{'August 2020'})
xline(157,'--k',{'Sept 2020'})
xline(187,'--k',{'Oct 2020'})
xline(217,'--k',{'Nov 2020'})
xline(247,'--k',{'Dec 2020'})
xline(277,'--k',{'Jan 2021'})
xticks([ 7 37 67 97 127 157 187 217 247 277])
xticklabels({'4/2020',' 5/2020','6/2020', '7/2020','8/2020','9/2020','10/2020','11/2020','12/2020','1/2021'})
yticks([ 3000 6000 9000 12000 15000 18000 21000])
yticklabels({'3000', '6000', '9000', '12000', '15000','18000','21000'})

figure(3)
plot(tdata, smooth_epiData,'o',epiForward,ppval(pp_smooth,epiForward))
hold on 
title({'B-Spline Interpolation of Florida Department of Health Data (With Smoothing)'}, 'Fontsize',20)
% hold on 
% plot(epiForward,ppval(p_der_smooth,epiForward),'-k')
hold on 
plot(tdata, epiData, 'r.', 'MarkerSize',20)
xlabel('Dates starting from March 23','fontweight','normal','fontsize',15)
ylabel('Number of Daily COVID-19 Incidences','fontweight','normal','fontsize',15)
xline(7,'--k',{'April 2020'})
xline(37,'--k',{'May 2020'})
xline(67,'--k',{'June 2020'})
xline(97,'--k',{'July 2020'})
xline(127,'--k',{'August 2020'})
xline(157,'--k',{'Sept 2020'})
xline(187,'--k',{'Oct 2020'})
xline(217,'--k',{'Nov 2020'})
xline(247,'--k',{'Dec 2020'})
xline(277,'--k',{'Jan 2021'})
xticks([ 7 37 67 97 127 157 187 217 247 277])
xticklabels({'4/2020',' 5/2020','6/2020', '7/2020','8/2020','9/2020','10/2020','11/2020','12/2020','1/2021'})
yticks([ 3000 6000 9000 12000 15000 18000 21000])
yticklabels({'3000', '6000', '9000', '12000', '15000','18000','21000'})


simulation_time = [epiForward  epiForward_vaccine];

results = [params(5)*y_r(:,2); params(5)*y_v(:,3)];


figure(9)
plot(simulation_time,results,'LineWidth',2)
hold on 
plot(tdata, epiData, 'r.', 'MarkerSize',20)
xlabel('Dates starting from March 23','fontweight','normal','fontsize',18)
xline(7,'--k',{'April 2020'})
xline(37,'--k',{'May 2020'})
xline(67,'--k',{'June 2020'})
xline(97,'--k',{'July 2020'})
xline(127,'--k',{'August 2020'})
xline(157,'--k',{'Sept 2020'})
xline(187,'--k',{'Oct 2020'})
xline(217,'--k',{'Nov 2020'})
xline(247,'--k',{'Dec 2020'})
xline(277,'--k',{'Jan 2021'})
xticks([ 7 37 67 97 127 157 187 217 247 277])
xticklabels({'Apr 2020',' May 2020','Jun 2020', 'Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020','Jan 2021'  })
yticks([ 3000 6000 9000 12000 15000 18000 21000])
yticklabels({'3000', '6000', '9000', '12000', '15000','18000','21000'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')

function dy = Model_SIR_ODE(t,y,k,epiForward,beta_t)
N = 21000000;
dy = zeros(4,1);

% beta = k(1);
gamma = k(2);
nu = k(3);
% betaE = k(4);
kk = k(5);
q =  k(4)/k(1);

beta = interp1(epiForward,beta_t,t);

S = y(1);
E = y(2);
I = y(3);
R = y(4);

dy(1) =  - beta ;
dy(2) =  beta - kk*E;
dy(3) = kk*E - (gamma + nu)*I;
dy(4) = gamma*I;

end

function dy = Model_SIR_ODE_vaccination(y,k,final_beta)
N = 21000000;

dy = zeros(5,1);

gamma = k(2);
nu = k(3);

kk = k(5);
q =  k(4)/k(1);

beta_v = final_beta/(0.0162*N);
beta_v = final_beta/(0.02*N);
%beta_v = final_beta/ 7.3526e+04;
vaccine_rate = 0.01;
vacc_efficacy = 0.9;

S = y(1);
V = y(2);
E = y(3);
I = y(4);
R = y(5);

dy(1) =  - beta_v*S.*(I./N + q*E./N) - vaccine_rate*S;
dy(2) =  vaccine_rate*S - (1-vacc_efficacy)*beta_v*V.*(I./N + q*E./N); 
dy(3) =  beta_v*S.*(I./N + q*E./N) + (1-vacc_efficacy)*beta_v*V.*(I./N + q*E./N) - kk*E;
dy(4) = kk*E - (gamma + nu)*I;
dy(5) = gamma*I;

end
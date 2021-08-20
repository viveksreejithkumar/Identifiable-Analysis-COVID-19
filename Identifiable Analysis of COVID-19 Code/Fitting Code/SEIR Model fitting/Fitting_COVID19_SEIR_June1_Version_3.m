clear all
close all
clc

global  epiData COVID_Deaths epiForward epiMeasure 

epiData = [631
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
3008
2054
2739
3266
3180
3584
3117
2617
1688
7440
2506
3602
3222
3705
2512
1880
1809
2330
2560
3691
3207
2418
1827];

COVID_Deaths = [30
40
59
38
49
41
34
32
37
44
43
41
46
32
36
41
22
37
39
38
34
41
36
36
25
39
40
27
37
35
30
41
38
26
31
32
27
24
33
27
30
37
25
33
29
45
32
30
24
38
36
51
58
40
54
45
52
50
47
53
68
74
69
83
74
99
104
123
114
106
121
131
146
128
146
165
158
150
164
179
161
156
163
160
180
163
149
155
158
167
164
165
151
164
181
183
154
125
126
130
147
142
133
134
114
112
89
94
91
99
91
90
73
71
62
72
78
53
65
52
44
51
51
40
38
39
31
39
36
24
24
17
6
2
1];

tdata = 0:1:134;

dt = 0.1; tf = 134;

epiForward = (0:dt:tf);
epiMeasure = 1:1/dt:length(epiForward);

params = [0.192190601961406,0.0952299449157147,0.00544027853949769,...
          1.41217937075999,0.197631293806486, 210000, 100, 10, 30];
%params = [0.0423722646592824,0.0178586841700095,0.0110252925203845,...
       %  1.787270389284732,0.00948158165896732, 210000, 100, 10, 2];
%params = [0.687747284192075,0.107133701152657,0.00611265128450433,0.000361870661769181,0.0448784080903576, 2100000, 1000, 100, 20];      

params = [0.319842880975307,0.101404098672993,0.00158432538080955,0.0311010086762042,0.0716314353709275];

params = [0.129384747432626,0.128305841202452,0.00161976378976192,0.142576732997668,0.144606587656209];

params = [0.0449279235259503,0.0895832792817166,0.00168484947105085,0.386591378082057,0.406022988314169];

params = [0.0238100043681300,0.0838460633128175,0.00168279866707977,0.510483671467554,0.489966077849274];

params = [0.0133183432504680,0.108953096839891,0.00172124565711584,0.557213670310253,0.499997101312672];

params = [0.0408828457254540,0.199999971959692,0.00212502382912483,0.507146513140100,0.499999990050984];

params = [0.0397588631240176,0.199999105214192,0.00212256973108893,0.505422687785417,0.499999994774246];

params = [0.0288252666791333,0.199999988628206,0.00206965866111908,0.524044990857231,0.499999992324709];

params = [0.0287508788221480,0.199451658704021,0.00207275032138021,0.527838280139914,0.499977658873889];

%lb = zeros(size(params));

lb = [0 1/30 0 0 1/15]
  
ub = [1e+5 0.2 1e+5 1e+5 0.5];

[params,fval] =  fminsearchbnd(@err_in_data,params,lb,ub,optimset('Display','iter'))

N = 650000*4.5;
initial_cond = [N-23 23 1300 32000];

[~, y_r] = ode15s(@(t,y)Model_SIR_ODE(y,params),epiForward,initial_cond);


 %w_incidence = 1/(mean(epiData)^2);
 %w_deaths = 1/(mean(COVID_Deaths)^2);
 Incidences = params(5)*y_r(epiMeasure(:),2);
 Deaths = params(3)*y_r(epiMeasure(:),3);
 
 error_in_incidences = ((Incidences - epiData)'*(Incidences - epiData))%*w_incidence;
 error_in_deaths = ((Deaths - COVID_Deaths)'*(Deaths - COVID_Deaths))%*w_deaths;   

figure(1)
plot(epiForward,params(5)*y_r(:,2),'LineWidth',2)
hold on 
plot(tdata, epiData, 'r.', 'MarkerSize',20)
title({'COVID-19 Incidences', 'May 2nd - September 13th'}, 'Fontsize',14)
xlabel('Days since May 2nd','Fontsize',14)
ylabel('Incidences','Fontsize',14)
set(gca,'linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',14)
xline(0,'--k',{'May 2nd, 2020'})
xline(20,'--k',{'May 22nd, 2020'})
xline(40,'--k',{'June 11th, 2020'})
xline(60,'--k',{'July 1st, 2020'})
xline(80,'--k',{'July 21st, 2020'})
xline(100,'--k',{'August 10th, 2020'})
xline(120,'--k',{'August 30th, 2020'})
xline(134,'--k',{'September 13th, 2020'})
xticks([0 20 40 60 80 100 120 134])
xticklabels({'5/2',' 5/22','6/11', '7/1','7/21','8/10','8/30','9/13'})

figure(2)
plot(epiForward,params(3)*y_r(:,3),'LineWidth',2.5)
hold on 
plot(tdata, COVID_Deaths, 'r.', 'MarkerSize',20)
title({'COVID-19 Deaths', 'May 2nd - September 13th'}, 'Fontsize',14)
xlabel('Days since May 2nd','Fontsize',14)
ylabel('Deaths','Fontsize',14)
set(gca,'linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',14)
xline(0,'--k',{'May 2nd, 2020'})
xline(20,'--k',{'May 22nd, 2020'})
xline(40,'--k',{'June 11th, 2020'})
xline(60,'--k',{'July 1st, 2020'})
xline(80,'--k',{'July 21st, 2020'})
xline(100,'--k',{'August 10th, 2020'})
xline(120,'--k',{'August 30th, 2020'})
xline(134,'--k',{'September 13th, 2020'})
xticks([0 20 40 60 80 100 120 134])
xticklabels({'5/2',' 5/22','6/11', '7/1','7/21','8/10','8/30','9/13'})


function error_in_data = err_in_data(k) 

global  epiData COVID_Deaths epiMeasure epiForward 

N = 650000*4.5;
initial_cond = [N-23 23 1300 32000];

 [~,y] = ode15s(@(t,y)Model_SIR_ODE(y,k),epiForward,initial_cond);
 
 Incidences = k(5)*y(epiMeasure(:),2);
 Deaths = k(3)*y(epiMeasure(:),3); 
 
 w_incidence = 1/((mean(epiData)^2)*length(epiData));
 w_deaths = 1/((mean(COVID_Deaths)^2)*length(COVID_Deaths));
%  
%  error_in_data = ((Incidences - epiData)'*(Incidences - epiData))*w_incidence+...
%                  ((Deaths - COVID_Deaths)'*(Deaths - COVID_Deaths))*w_deaths;           
  error_in_data = sum((Incidences - epiData).^2)*w_incidence+...
                 sum((Deaths - COVID_Deaths).^2)*w_deaths;   
end

function dy = Model_SIR_ODE(y,k)
N = 650000*4.5;

dy = zeros(4,1);

beta = k(1);
gamma = k(2);
nu = k(3);
betaE = k(4);
k = k(5);

S = y(1);
E = y(2);
I = y(3);
R = y(4);

dy(1) =  - betaE*S.*E./N - beta*S.*I./N;
dy(2) = beta*(S).*I./N + betaE*(S).*E./N - k*E;
dy(3) = k*E - (gamma + nu)*I;
dy(4) = gamma*I; 
end
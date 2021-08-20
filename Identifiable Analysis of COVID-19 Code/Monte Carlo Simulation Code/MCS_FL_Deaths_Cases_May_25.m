clear all
close all
clc

global epiForward initial_cond epiMeasure 

numiter = 1000;

true_params = [0.192191938697128,0.265372537962109,0.456404771410813,0.0952279292130327,0.00544021233585557,0.00587553500031035,1.41217214680023,0.00828729549913443,0.197630685802428];


dt = 0.1; tf = 61;

epiForward = (0:dt:tf);
epiMeasure = 1:1/dt:length(epiForward);

N = 210000;
initial_cond = [N-193 0 193 0 0];

[~, y_trp] = ode15s(@(t,y)Model_SIR_ODE(y,true_params),epiForward,initial_cond);

 Incidences_trp = true_params(9)*y_trp(epiMeasure(:),3);
               
 Deaths_trp = true_params(5)*y_trp(epiMeasure(:),4); 


X = zeros(9,numiter);


noiselevel = [0, 0.01, 0.05, 0.1, 0.2];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:5
    
rng default
noiselev = noiselevel(noisei)

for i = 1:numiter
    i
    
  epiData =  (noiselev*(Incidences_trp).*randn(length(epiMeasure),1)) + Incidences_trp;
  deathData =  (noiselev*(Deaths_trp).*randn(length(epiMeasure),1)) + Deaths_trp;

   k = true_params;
   lb = zeros(size(k));
   ub = [1e+4 1 1 0.2 1 1e+3 1e+4 1 0.5];
   
   [k,~] =  fminsearchbnd(@(k)err_in_data(k,epiData,deathData),k,lb,ub,optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-8,'TolX',1e-8))
 
   X(:,i) = k';
end

arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
        
    end
    
    total_ARE(noisei,:) = arescore;
end

save('MCS_May_25')
function error_in_data = err_in_data(k,epiData,deathData) 

global epiForward initial_cond epiMeasure 

 [~,y] = ode15s(@(t,y)Model_SIR_ODE(y,k),epiForward,initial_cond);
 

 Model_Predictions_epi = k(9)*y(epiMeasure(:),3);
 Model_Predictions_death = k(5)*y(epiMeasure(:),4); 
 
 error_in_data = (Model_Predictions_epi - epiData)'*(Model_Predictions_epi - epiData) +...
     ((Model_Predictions_death - deathData)'*(Model_Predictions_death - deathData))*100;      
  
end

function dy = Model_SIR_ODE(y,k)
N = 210000;


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
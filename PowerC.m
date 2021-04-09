function [AEPtotal,CF,LCOE,WTPerf,OpRange,KeyMet,Farm,U_in,FarmPerf] = PowerC(A,k,type)
%4: Constrained power curve function based on Betz efficiency from BEM
%solver
%Inputs of wind turbne flow control configuration and turbine wind
%distribution data
%Outputs:
    %Idealised Power curves
    %ind turbine performance curves
    %AEP, CF & LCOE
WindTurbineModeller 

    global rho R r
    
    MinV = 0.5; %Min windspeed
    MaxV = 80; %Max windspeed
    MinU = 3; %Min windspeed (farm calcs)
    MaxU = 35.5; % Max windspeed (farm calcs)
    P_in = 42900; %Power at cut in (Clean config)
    P_rated = 5.2966e6; %Actual Power to give rated power at Generator (losses)
    P_r = 5e6; %Rated 
    P_cout = 63.22e6; %Power at cut-out on idealised curve (Clean config)
    
    %Velocity ranges
    du = 0.5;
    dv = 1;
v1 = MinV + dv/2;
v2 = MaxV - dv/2;
u1 = MinU + du/2;
u2 = MaxU - du/2;
n = 1;
s = 1;
t = 1;
%Rotor Area
Area = pi*(R^2-r(1)^2);

%Betz Efficiency Calculation
[ ~, Eff, ~, ~] = VRange(A, k, 1.267, 3, 30, type);

for V = [v1:dv:v2]
%idealised Power Calculation
    P_i = (Eff/100) * 0.5 * rho * Area * V^3;
    IdealisedP(s,:) = P_i;
    Windspeed(s,:) = V;
    
    s=s+1;
end

%Interpolate Idealised Power to calculate operating velocities and display
 
V_in = interp1(IdealisedP,Windspeed,P_in,'pchip');
V_r = interp1(IdealisedP,Windspeed,P_rated,'pchip');    
V_cout = interp1(IdealisedP,Windspeed,P_cout,'pchip');

OpRange = [V_in V_r V_cout];

%Manual overide of velocities to calculated combined solutions

%FF combined Op range
% V_in = 2;
% V_r = 10.0; 
% V_cout = 71.3;
%VGs & FFs combined Op range
% V_in = 1.9;
% V_r = 9.5; 
% V_cout = 49.4;

%FF Location Op range
% V_in = 2.20;
% V_r = 10.94;
% V_cout = 40.27;

%FF 40-60% & VGs Op range
% V_in = 2.03219026370313;
% V_r = 10.1194215844166;
% V_cout = 27.9851;

for V = [v1:dv:v2]
    
    
    
    %Probability
    Prob = exp(-((V-dv/2)/A)^k)-exp(-((V+dv/2)/A)^k);
    
     %Normalised Power P_n
   q = 2; %Power
   P_n = IdealisedP(n)*(V^q - V_in^q)/(V_r^q - V_in^q);
   
   %Create Power curve that is limited based on cut off speeds (baseline)
   
   if V < V_in
       P_e =0;
   elseif V >= V_in & V < V_r
       P_e = P_n;
   elseif V >= V_r & V<= V_cout
       P_e = P_r;
   elseif V > V_cout
       P_e = 0;
   end
   
   %Annual Energy Production Calculation
   AEP = P_e * Prob * 8760;
   
   WTPerf(n,:) = [V P_i P_e Prob AEP];
  
   n = n+1;
end

% Wind Farm Performance calculations (4 x 4 wind farm)
for U = [u1:du:u2]
    % Calculates the velocities at each turbine row (wake calculations
    [~, ~, ~, U_w,CT] = SingleVelocity(U, 1.267, 3, type);
    
    for i= 1:length(U_w) %loop through 4 velocities
        
        Prob_f = exp(-((U-du/2)/A)^k)-exp(-((U+du/2)/A)^k);
        P_i_f = (Eff/100) * 0.5 * rho * Area * U_w(i)^3;
        P_n_f = P_i_f*(U_w(i)^q - V_in^q)/(V_r^q - V_in^q);
        
        if U_w(i) < V_in
            P_e_f =0;
        elseif U_w(i) >= V_in & U_w(i) < V_r
            P_e_f = P_n_f;
        elseif U_w(i) > V_cout | U_w(1) > V_cout
            P_e_f = 0;
        elseif U_w(i) >= V_r & U_w(i)<= V_cout 
            P_e_f = P_r;
        
        end
        Farm_p(:,i) = P_e_f;
        AEP_f = P_e_f * Prob_f * 8760;
        AEP_Farm(:,i) = AEP_f;
        
    end
    Farm(t,:) = AEP_Farm;
    FarmPower(t,:) = Farm_p;
    U_out(t,:) = U;
    U_in(t,:) = U_w;
    CTmat(t,:) = CT;
    t = t+1;
    
end
   %Annual Energy Production
   AEPtotal = sum(WTPerf(:,5))/1e6
   
   %Annual Energy Production (Farm)
   AEPFarmTotal = (4*sum(sum(Farm)))/1e6;
   
   %Farm Efficiency
   FarmEff = AEPFarmTotal/(16*sum(Farm(:,1))/1e6);
   
   %Capacity Factor Calculation
   CF = sum(WTPerf(:,5))/(P_r*8760);
   
   %Farm Capacity factor
   FarmCF = AEPFarmTotal/((16*P_r*8760)/1e6);
   
   %sLCOE Calculation simple Levelized Cost of Electricity
   
   i = 0.1; %Interest rate (percent)
   t = 20; %number of years
   CRF = (i*(1+i)^t)/((1+i)^t - 1); %Capital Recovery Factor
   C_o = 2950; %Overnight Capital Cost ($/kW)
   O = 95; %Fixed operating cost ($/kW-yr)
   
   LCOE = (C_o*CRF+O)/(8760*CF) %($/kWh)
   FarmLCOE = (C_o*CRF+O)/(8760*FarmCF);
   
   %figure(1)
   %plot(WTPerf(:,1),WTPerf(:,3)); %Performance Curve
    %figure(2)
%     plot(WTPerf(:,1),WTPerf(:,5)); %AEP
   %plot(WTPerf(:,1),WTPerf(:,2));
   %plot(Windspeed,IdealisedP)
   KeyMet = [Eff AEPtotal CF LCOE V_in V_r V_cout]; %Key output metrics
   FarmPerf = [AEPFarmTotal FarmEff FarmCF FarmLCOE]; 
end


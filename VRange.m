function [Out, Eff, V_out,FarmVs] = VRange(A, k, omega, MinV0, MaxV0, type)
%3: ANNUAL ENERGY - loop WTSingleVelocity to find the moments across the
%entire velocity range. Combine this with the frequency information to get
%the AEP.

global rho R r

R = 64.4;
B = 3;

dv = 1;
v1 = MinV0 + dv/2;
v2 = MaxV0 - dv/2;
n = 1;

for V0 = [v1:dv:v2]
   %Moments
   
   [MT, MN, ~, U_w] = SingleVelocity(V0, omega, B, type);
   
   
   %Power calc
   
       P = sum(MT) * B * omega;
       if P < 0 
           P=0;
       end
   
   %Torque Calc
   
   Q = (sum(MT)*B)/1e6;
   
   %Ideal Power
   
   Area = pi*(R^2-r(1)^2);
   
   P_ideal = 16/27 * 0.5 * rho * V0^3 * Area;
   
   %Probability
   Prob = exp(-((V0-dv/2)/A)^k)-exp(-((V0+dv/2)/A)^k);
   
   %AEP
   AEP = P * Prob * 8760;
   
   %AEP Ideal
   AEP_ideal = P_ideal * Prob * 8760;
   
   
   %Output
   Out(n,:) = [ P P_ideal Prob AEP AEP_ideal ]; 
   FarmVs(n,:) = [U_w];
   V_out(n,:) = V0;
  
   n = n+1;
end


AEP = sum(Out(:,4));
AEP_ideal = sum(Out(:,5));

%Turbine Efficiency
Eff = (AEP/AEP_ideal)*100


% plot(V_out,Out(:,1));
%figure
%plot(V_out,PowerC(:,1));

end
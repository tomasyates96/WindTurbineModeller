function [MT, MN, Result, U_w,CT] = SingleVelocity(V0, omega,B, type)
%2: WHOLE ROTOR - loop WTInducedCalcs to find the values for all radii,
%then integrate these to get the normal and tangential moment at the blade
%root.

n = 1;

global rho mu R r dr Twist Chord 

for i = 1:length(r);
    
    %Initial guesses
    a = 0;
    adash = 0;
    
    %Theta calculation
    
    Theta = Twist(n)*(pi/180);
   
    %equation for Chord
    C = Chord(n);
    
    %node position
    y = r(n);
    
    %node size
    Ns = dr(n);
    
    [a, adash, phi, Cn, Ct, Cl, Cd, Re] = InductionFactors(a, adash, V0, omega, Theta, C, B, y, type);
    
    Result(n,:) = [y, a, adash, phi, Cn, Ct, Cl, Cd, Re];
    
    %Vrel
    Vrel = ((V0*(1-a))^2 + (omega*r(n)*(1+adash))^2)^0.5;
    
    %Equation for MT
    
    MT(n,:) = (0.5*rho*Vrel^2*C*Ct)*Ns*y;
    
    %Equation for MN
    
    MN(n,:) = (0.5*rho*Vrel^2*C*Cn)*Ns*y;
    
    
    n = n+1;
end

%Wake Velocity calculations for downwind turbines
CT = mean(Result(:,6));
k = 0.05; %wake decay coefficient
x = 10*(2*R); %turbine seperation distance (10D)

U_w1 = V0; %Inflow velocity at row 1

U_w12 = U_w1*(1-(1-sqrt(1-CT))*(R/(R+k*x))^2); %Wake 1 - 2 (turbine 2 inflow)

U_w13 = U_w1*(1-(1-sqrt(1-CT))*(R/(R+k*2*x))^2); % 1 - 3

U_w23 = U_w12*(1-(1-sqrt(1-CT))*(R/(R+k*x))^2); % 2 - 3

U_w3 = U_w1*(1-sqrt((1-U_w13/U_w1)^2 + (1-U_w23/U_w1)^2)); %(turbine 3 inflow)

U_w14 = U_w1*(1-(1-sqrt(1-CT))*(R/(R+k*3*x))^2); % 1 - 4

U_w24 = U_w12*(1-(1-sqrt(1-CT))*(R/(R+k*2*x))^2); % 2 - 4

U_w34 = U_w3*(1-(1-sqrt(1-CT))*(R/(R+k*x))^2); % 3 - 4

U_w4 = U_w1*(1-sqrt((1-U_w14/U_w1)^2 + (1-U_w24/U_w1)^2 + (1-U_w34/U_w1)^2)); %(turbine 4 inflow)

U_w = [ U_w1 U_w12 U_w3 U_w4 ];

end
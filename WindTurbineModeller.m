
%Script to run all calculation functions to model NREL wind turbine
%Constants

global mu rho Chord Twist r B dr R VGalpha VGcl FFPos_alpha FFPos_Cl FFPos_Cd_alpha ...
    FFPos_Cd FFNeg_alpha FFNeg_Cl FFNeg_Cd_alpha FFNeg_Cd

mu = 1.81e-5; %Dynamic viscosity of air
rho = 1.225; % Air Density at SL
B = 3; %Number of Blades
R = 63; %Blade Radius
load('FYP_NREL_5MW_Input_Data.mat');
load('VGdata.mat');
load('FFPos_Cl');
load('FFPos_Cd');
load('FFNeg_Cl');
load('FFNeg_Cd');

%-----------
% Chord Data
%-----------
Chord = FYPdata(:,5);

%-----------
% Twist Data
%-----------
Twist = FYPdata(:,3);

%-----------
% Radial Position Data
%-----------
r = FYPdata(:,2);

%----------
%Node Size
%----------
dr = FYPdata(:,4);

%----------
%VG data
%----------
VGalpha = FYPdata1(:,1);
VGcl = FYPdata1(:,2);

%----------
%FF data
%----------

%Max posritive deflection - lift increasing
FFPos_alpha = FYPdata2(:,2);
FFPos_Cl = FYPdata2(:,3);

FFPos_Cd_alpha = FYPdata3(:,2);
FFPos_Cd = FYPdata3(:,3);

%Max negative deflection - lift decreasing
FFNeg_alpha = FYPdata4(:,2);
FFNeg_Cl = FYPdata4(:,3);

FFNeg_Cd_alpha = FYPdata5(:,2);
FFNeg_Cd = FYPdata5(:,3);








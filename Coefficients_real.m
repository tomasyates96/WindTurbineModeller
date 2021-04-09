function [Cl, Cd] = Coefficients_real(Alpha, Re) %Add in airfoil to input
%CALCULATE Cl and Cd. 
persistent FL FD
%For alpha < 20 this is done by an approximation to the experimental
%data. The function TriScatteredInterp fits a surface to the raw data which
%then allows values to be extracted. The more data, the more accurate the
% method.

%For alpha > 20deg it is an algebraic method based on the methods described
%in the notes and Leishman.

Alpha=Alpha*180/pi; %Convert from radians into degrees.
if Re < 750000 %Lower Limit of Experimental Data
    Re = 750000;
 elseif Re> 1500000 %Upper Limit of Experimental Data
    Re = 1500000;
end
if Alpha<=40 && Alpha>=-20
    if isempty(FL) %Check if the fit already exists to avoid repetition
        %IMPORT DATA
            
        LiftData=dlmread('S809_Lift_Real.txt');
        DragData=dlmread('S809_Drag_Real.txt');
        
      
        %FIT SURFACE
        FL = TriScatteredInterp(LiftData(:,1),LiftData(:,2),LiftData(:,3)); 
        FD = TriScatteredInterp(DragData(:,1),DragData(:,2),DragData(:,3));
        
        PlotFit=0;
        if PlotFit == 1
            %PLOT LIFT
            [qx,qy] = meshgrid(10000: 100000: 6000000, -20:1:20);
            qz = FL(qx,qy);
            mesh(qx,qy,qz);
            hold on;
            plot3(LiftData(:,1),LiftData(:,2),LiftData(:,3),'ko');
            plot3(Re, Alpha, FL(Re, Alpha), 'ro', 'MarkerFaceColor', 'r') 
            
            %PLOT DRAG
            figure
            qz = FD(qx,qy);
            mesh(qx,qy,qz);
            hold on;
            plot3(DragData(:,1),DragData(:,2),DragData(:,3),'ko');
            plot3(Re, Alpha, FD(Re, Alpha), 'ro', 'MarkerFaceColor', 'r') 
        end
    end
    
    Cl=FL(Re, Alpha);
    Cd=FD(Re, Alpha);
    
    if Alpha > 35 || Alpha < -15
       %Blend with highAlpha result to avoid discontinuities at alpha=20. 
       %Discontinuities can cause instability in iterative solutions.
       l0 = 0.1867;
       l1 = 1.4885;
       l2 = 0.1991;
       Cl=(Cl*(20-abs(Alpha))/5)+((l0+l1*sin(2*Alpha*pi/180)+l2*sin(4*Alpha*pi/180))*(abs(Alpha)-15)/5);
       
       d0 = 1.1657;
       d1 = -1.0058;
       d2 = -0.1253;
       Cd=(Cd*(20-abs(Alpha))/5)+((d0+d1*cos(2*Alpha*pi/180)+d2*cos(4*Alpha*pi/180))*(abs(Alpha)-15)/5);
    end
end

%USE High Alpha Solution if alpha is greater than 20 or if the
%approximation method failed. Contants from Apostolyuk
if Alpha > 40 || Alpha < -20 || exist('Cd', 'var') == 0 ||isnan(Cd) == 1
    d0 = 1.1657;
    d1 = -1.0058;
    d2 = -0.1253;
    Cd=d0+d1*cos(2*Alpha*pi/180)+d2*cos(4*Alpha*pi/180);
end

if Alpha > 40 || Alpha < -20 || exist('Cl', 'var') == 0 || isnan(Cl) == 1
    l0 = 0.1867;
    l1 = 1.4885;
    l2 = 0.1991;
    Cl=l0+l1*sin(2*Alpha*pi/180)+l2*sin(4*Alpha*pi/180);
end
end

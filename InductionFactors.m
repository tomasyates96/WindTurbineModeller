function [a, adash, phi, Cn, Ct, Cl, Cd, Re] = InductionFactors(a, adash, V0, omega, theta, Chord, B, y, type)
%1: SINGLE ELEMENT: use an iterative solution to find the values of a,
%adash, phi, Cn and Ct at a particular radius.

%Constants
global rho mu R r VGalpha VGcl FFPos_alpha FFPos_Cl FFPos_Cd_alpha ...
    FFPos_Cd FFNeg_alpha FFNeg_Cl FFNeg_Cd_alpha FFNeg_Cd



% Relaxation factor 
k = 0.1;

%a_c for Glauert correction (Spera 1994)
a_c = 0.2;

Counter = 0;
Error = 1;
tol = 0.001;
    
    
    
    while Error > tol & Counter < 200 
        
        n=1;
        %Flow angle
        tan_phi = ((1-a)*V0)/((1+adash)*omega*y);
        phi = atan(tan_phi);
        
        %AoA
        alpha = phi - theta;
       
        
        %Relative velocity
        Vrel = ((V0*(1-a))^2 + (omega*y*(1+adash))^2)^0.5;
        
        
        if ~isreal(alpha)
            disp(alpha);
        end
        
        %Reynold's number
        Re = (rho*Vrel*Chord)/mu;
        
        %Coefficients of lift and drag
        [Cl, Cd] = Coefficients_real(alpha,Re);
 
        
        %Conversion of force coefficients
        Cn = (Cl)*cos(phi) + Cd*sin(phi);
        Ct = ((Cl)*sin(phi) - Cd*cos(phi));
        
        %Solidity
        
        sigma = (B*Chord)/(2*pi*y);
        
        %Prandtl's Tip Loss Correction Factor
        
        f = (B/2)*((R-y)/(y*sin(phi)));
        F = (2/pi)*acos(exp(-f));
        if  F < 0 || F > 1 || ~isreal(F)
            F = 1;
        end
        
        % a & adash calculation with Prandtl's tip loss
        %Glauert Correction also included for a values above a_c = 0.2
        
        K = (4*F*sin(phi)^2)/(sigma*Cn);
        
        a_out = 1/(((4*F*sin(phi)^2)/(sigma*Cn))+1);
        
            if a_out > a_c
                % for higher values of a, Glauert correction applied
                a_out = 0.5 *(2+ (K*0.6)-sqrt((K*0.6+2)^2+4*(K*0.2^2 -1)));
            end
            
        
        %a_out =real(a_out);
        
        if ~isreal(a_out)
            disp(a_out)
        end
        
        %a dash calculation
        if Counter >100
            adash_out = 0;
            
        else
            adash_out = 1/(4*sin(phi)*cos(phi)/(sigma*Ct)-1);
        end
        
        Error = abs(a_out-a) + abs(adash_out-adash);
        
        %relaxation factors
        a_out = k*(a_out-a) + a;
        adash_out = k*(adash_out-adash) + adash;
        
        a = a_out;
        
        if Counter > 100
            adash = 0;
        else
        adash = adash_out;
        end
        
        Counter = Counter +1;
        
        
    end
    
     %---------------------------------
        % Lift Coefficient Modification
     %---------------------------------
     if type == 'A' %A is clean configuration (no mods)
         X=0;
         D=0;
         Y=0;
     %---------  Vortex Generators ----------------
     elseif type == 'B' %B is VG generators only
         D=0;
         Y=0;
         if alpha > 0.069661326 && alpha < 0.38503185  % Limits of data available
             X = interp1(VGalpha,VGcl,alpha,'linear'); %Interpolate VG data for Cl delta
         else
             X=0;
         end
     %---------  FF (max positive deflection) ----------------
     elseif type == 'C' %&& y > 0.8*R
         Y=0;
             if alpha > -0.162347490035133 && alpha < 0.394433162856721 % Limits of data available
                 X = interp1(FFPos_alpha,FFPos_Cl,alpha,'linear'); %Interpolate FF data for Cl delta
             else
                 X=0;
             end
             
             if  alpha > -0.163472299434915 && alpha < 0.214838595358271
                 D = interp1(FFPos_Cd_alpha,FFPos_Cd,alpha,'linear'); %Interpolate FF data for Cd delta
             else
                 D =0;
             end

     %---------  FF (max negative deflection) ---------------- 
     elseif type == 'D' %&& y > 0*R && y < 0.2*R%D is FF (max negative deflection)
         Y=0;
         if alpha > -0.163054496756439 && alpha < 0.288861113134950 % Limits of data available
             X = interp1(FFNeg_alpha,FFNeg_Cl,alpha,'linear'); %Interpolate FF data for Cl delta
         else
             X=0;
         end
         
         if  alpha > -0.163054496756439 && alpha < 0.210944628670917  % Limits of data available
             D = interp1(FFNeg_Cd_alpha,FFNeg_Cd,alpha,'linear'); %Interpolate FF data for Cd delta
         else
             D =0;
         end
     %---------  FF & VGs (max positve deflection) ----------------     
     elseif type == 'E' %E is VGs and FF (max positive)
         
         if alpha > 0.069661326 && alpha < 0.38503185  % Limits of data available
             X = interp1(VGalpha,VGcl,alpha,'linear'); %Interpolate VG data for Cl delta
         else
             X=0;
         end
         if alpha > -0.162347490035133 && alpha < 0.394433162856721 %&& y > 0.4*R && y < 0.6*R % Limits of data available
             Y = interp1(FFPos_alpha,FFPos_Cl,alpha,'linear'); %Interpolate FF data for Cl delta
         else
             Y=0;
         end
         
         if  alpha > -0.163472299434915 && alpha < 0.214838595358271 %&& y > 0.4*R && y < 0.6*R
             D = interp1(FFPos_Cd_alpha,FFPos_Cd,alpha,'linear'); %Interpolate FF data for Cd delta
         else
             D =0;
         end
     
     %---------  FF & VGs (max negative deflection) ---------------- 
     elseif type == 'F' %F is VGs and FF (max negative)
         if alpha > 0.069661326 && alpha < 0.38503185  % Limits of data available
             X = interp1(VGalpha,VGcl,alpha,'linear'); %Interpolate VG data for Cl delta
         else
             X=0;
         end
         if alpha > -0.163054496756439 && alpha < 0.288861113134950 %&& y > 0.4*R && y < 0.6*R % Limits of data available
             Y = interp1(FFNeg_alpha,FFNeg_Cl,alpha,'linear'); %Interpolate FF data for Cl delta
         else
             Y=0;
         end
         
         if  alpha > -0.163054496756439 && alpha < 0.210944628670917 %&& y > 0.4*R && y < 0.6*R  % Limits of data available
             D = interp1(FFNeg_Cd_alpha,FFNeg_Cd,alpha,'linear'); %Interpolate FF data for Cd delta
         else
             D =0;
         end
         
     else
         X=0;
         Y=0;
         D=0;
     end
    Cn = (Cl+X+Y)*cos(phi) + (Cd+D)*sin(phi);
    Ct = ((Cl+X+Y)*sin(phi) - (Cd+D)*cos(phi));
end


function [msrs_out] = Nom_Msrs(C, t, i)
    %==========================================================================
    % [msrs_out] = Nom_msrs(C, t, i)
    %
    %  Finds the nominal measurements at a given time for a given station
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for   int [1-12] 
    %
    % OUTPUT:       
    %    
    %                  - Range from station to nominal trajectory     km
    %                                     
    % Coupling:
    %
    %  X_i                - Gives the X location of station i at a given
    %                       time
    %  Y_i                - Gives the Y location of station i at a given
    %                       time
    %
    %==========================================================================
    msrs_out = [Rho_Nom(C, t, i); ((X_Nom(C, t) - X_i(C, t, i)) * (X_Nom_Dot(C, t) - dX_i(C, t, i)) + ...
                (Y_Nom(C, t) - Y_i(C, t, i)) * (Y_Nom_Dot(C, t) - dY_i(C, t, i))) / Rho_Nom(C, t, i); 
                atan2((Y_Nom(C, t) - Y_i(C, t, i)),(X_Nom(C, t) - X_i(C, t, i)))];

end


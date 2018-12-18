 function phi = Phi_I(C, t, i) 
    %==========================================================================
    % [phi] = Phi_I(C, t, i)
    %
    % Finds the nominal elevation angle of the spacecraft relative to a
    % given ground station 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for   int [1-12] 
    %
    % OUTPUT:       
    %    
    %  phi               - Nominal elevation angle of the spacecraft     rad 
    %                                     
    % Coupling:
    %
    %  Y_i               - Finds station y position
    %  X_i               - Finds station x position
    %
    %========================================================================== 
    phi = atan2((Y_Nom(C, t) - Y_i(C, t, i)), (X_Nom(C, t) - X_i(C, t, i)));
 end

function H_out = H_i(C, t, i)
    %==========================================================================
    % [H_out] = H_i(C, t, i)
    %
    % Consructs the observation matrix for a single station at a given time
    % 
    % INPUT:               Description                                   Units
    %
    %  C                  - Constants Object                              NA
    %  time               - time at which to initialize h matrix          sec
    %  i                  - number of station to generate h matrix for    int [1-12] 
    %
    % OUTPUT:       
    %    
    %  H_out               - Observation matrix for station i             (nx25)
    %                                     
    % Coupling:
    % 
    %  rho_nom             - evaluates the nominal range from station i   
    %  Q_nom               - evaluates a nasty math expression I
    %                        arbitrarily called Q
    %
    %==========================================================================    
    H_out = [(X_Nom(C, t) - X_i(C, t, i)) / Rho_Nom(C, t, i), 0, ... 
             (Y_Nom(C, t) - Y_i(C, t, i)) / Rho_Nom(C, t, i), 0;   
             (Y_i(C, t, i) - Y_Nom(C, t)) * Q_Nom(C, t, i) / Rho_Nom(C, t, i)^3, ...  
             (X_Nom(C, t) - X_i(C, t, i))/ Rho_Nom(C, t, i), ...
             (X_Nom(C, t) - X_i(C, t, i)) * Q_Nom(C, t, i) / Rho_Nom(C, t, i)^3, ...  
             (Y_Nom(C, t) - Y_i(C, t, i)) / Rho_Nom(C, t, i);        
            -(Y_Nom(C, t) - Y_i(C, t, i)) / Rho_Nom(C, t, i)^2, 0, ...
             (X_Nom(C, t) - X_i(C, t, i)) / Rho_Nom(C, t, i)^2, 0];
end    


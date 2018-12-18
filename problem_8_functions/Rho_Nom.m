function rho = Rho_Nom(C, t, i) 
    %==========================================================================
    % [rho] = Rho_Nom(C, t, i)
    %
    %  Finds the nominal range between the station i and satellite at a
    %  given time 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for   int [1-12] 
    %
    % OUTPUT:       
    %    
    %  rho                - Range from station to nominal trajectory     km
    %                                     
    % Coupling:
    %
    %  X_i                - Gives the X location of station i at a given
    %                       time
    %  Y_i                - Gives the Y location of station i at a given
    %                       time
    %
    %==========================================================================
    rho = sqrt((X_Nom(C, t) - X_i(C, t, i))^2 + (Y_Nom(C, t) - Y_i(C, t, i))^2);
end


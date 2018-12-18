function q = Q_Nom(C, t, i) 
    %==========================================================================
    % [q] = Q_nom(C, time)
    %
    % Solves a gnarly bit of math necessary to generate the H matrix
    %
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for    int [1-12] 
    %
    % OUTPUT:       
    %    
    % q                 - some resulting value from knarly math
    %                                     
    % Coupling:
    %
    % X_i            - Finds the x position of a station 
    % Y_i            - Finds the y position of a station
    % dX_i           - Finds the x velocity of a station
    % dY_i           - Finds the y velocity of a station
    %
    %========================================================================== 
    q = (X_i(C, t, i) - C.r0 * cos(C.n * t)) * (dY_i(C, t, i) - C.n * C.r0 * cos(C.n * t)) - ...
        (Y_i(C, t, i) - C.r0 * sin(C.n * t)) * (dX_i(C, t, i) + C.n * C.r0 * sin(C.n * t));
end

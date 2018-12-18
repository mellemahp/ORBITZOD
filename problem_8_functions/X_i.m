function [x_i_out] = X_i(C, t, i) 
    %==========================================================================
    % [x_i_out] = X_i(C, t, i)
    %
    % Gives the inertial X position of a given station 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for   int [1-12] 
    %
    % OUTPUT:       
    %    
    %  x_i_out           - x position of station in inertial frame       km
    %                                     
    % Coupling:
    %
    %  None 
    %
    %========================================================================== 
    x_i_out = C.Re * cos(C.omega_e * t + C.station_angles_0(i));
end

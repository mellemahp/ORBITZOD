function [y_i_out] = Y_i(C, t, i)
    %==========================================================================
    % [y_i_out] = Y_i(C, t, i)
    %
    % Gives the inertial y position of a given station 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to evaluate position          sec
    %  i                  - number of station to generate value          int [1-12] 
    %
    % OUTPUT:       
    %    
    %  y_i_out           - y position of station in inertial frame       km
    %                                     
    % Coupling:
    %
    %  None 
    %
    %==========================================================================    
    [y_i_out] = C.Re * sin(C.omega_e * t + C.station_angles_0(i));
end
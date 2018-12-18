function [dx_i_out] = dX_i(C, t, i)
    %==========================================================================
    % [dx_i_out] = dX_i(C, t, i)
    %
    % Gives the inertial X veloity of a given station 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to evaluate function          sec
    %  i                  - number of station to generate value          int [1-12] 
    %
    % OUTPUT:       
    %    
    %  x_i_out           - x velocity of station in inertial frame       km/sec
    %                                     
    % Coupling:
    %
    %  None 
    %
    %========================================================================== 
    dx_i_out = -C.Re * C.omega_e * sin(C.omega_e * t + C.station_angles_0(i));
end


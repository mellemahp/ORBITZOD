function [dy_i_out] = dY_i(C, t, i) 
    %==========================================================================
    % [dy_i_out] = dY_i(C, t, i)
    %
    % Gives the inertial y velocity of a given station 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to evaluate position          sec
    %  i                  - number of station to generate value for      int [1-12] 
    %
    % OUTPUT:       
    %    
    %  dy_i_out           - y velocity of station in inertial frame      km/sec
    %                                     
    % Coupling:
    %
    %  None 
    %
    %==========================================================================  
    dy_i_out = C.Re * C.omega_e * cos(C.omega_e * t + C.station_angles_0(i));
end

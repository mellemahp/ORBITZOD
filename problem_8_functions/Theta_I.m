function [stn_theta] = Theta_I(C, t, i) 
    %==========================================================================
    % [stn_theta] = Theta_I(C, time)
    %
    % Gives the stations angle in the inertial frame 
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for   int [1-12] 
    %
    % OUTPUT:       
    %    
    %  stn_theta      - angle off the X axis at time                    rad
    % 
    %                                     
    % Coupling:
    %
    %   X_i  -provides station's X position in inertial frame
    %   Y_i  -provides station's Y position in inertial frame  
    %
    %========================================================================== 
    stn_theta = atan2(Y_i(C, t, i) , X_i(C, t, i));
end


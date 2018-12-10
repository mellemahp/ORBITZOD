
clc 
clear all
close all 

%% Constants

C.mu = 3.986e5; % km^3/sec^2
C.r0 = 6678; % km 
C.period = 2 * pi * sqrt(C.r0^3 / C.mu); % sec
C.n = 2 * pi / C.period;
C.delta_t = 10; % sec 
C.Re = 6378; %km 
C.omega_e = 2 * pi / 86400; % earth rotation rate
C.station_angles_0 = arrayfun(@(i) (i - 1) * pi / 6, 1:12);
                                                   
%% Generate Truth States

% create the time vector for two periods of the orbit
times = 1:C.delta_t:14000; %C.period * 2;

% initialize the initial conditions of the orbit and perturbation vector
state = [C.r0; 0; 0; C.r0 * sqrt(C.mu / C.r0^3)];
p_vec= [0., 0.075, 0, -0.021]';

% set up nominal trajectory 
x_nom = C.r0 * cos(C.n * times);
y_nom = C.r0 * sin(C.n * times); 

% Run the linearized dynamics for the perturbations
for k = 1:(length(times)-1)
    p_vec(:, k+1) = F_tilde(C, times(k)) * p_vec(:, k);
end 

%% Plotting the linearized dynamics 
lin_dyn_titles = ["X", "Xdot", "Y", "Ydot"];

for k = 1:4
   subplot(4,1,k)
   plot(times, p_vec(k, :))
   title(strcat(lin_dyn_titles(k), " Perturbations"))
   if k == 2
       ylabel("Perturbation Magnitude (km)               ")
   end 
end

xlabel("Time (seconds)")

%% Generate Simulated Ground Station Data!

msrs = [];
for t_index = 1:(length(times)) 
    h_mat = Build_Full_H_Matrix(C, times(t_index));
    msrs(:, t_index) = h_mat * p_vec(:, t_index);
end

%% Plot Measurements

% F
x = 1;
for i = 1:12
    figure()
    subplot(3, 1, 1)
    plot(times, msrs(x,:))
    title(strcat('Station', num2str(i)))
    ylabel('Range (km)')
    subplot(3, 1, 2) 
    plot(times, msrs(x + 1,:))
    ylabel('Range Rate (km/s)')
    subplot(3, 1, 3)
    plot(times, msrs(x + 2,:))
    ylabel('Phi (rad)')
    xlabel('Time (sec)')
    x = x + 3;
end 


                    
%% FUNCTIONS
% ====================================================================================


%% Dynamics Matrix 

function [f_tilde_out] = F_tilde(C, t)
    %==========================================================================
    % [f_tilde_out] = F_tilde(C, time)
    %
    % Consructs the linearized dynamics matrix (linearized about the
    % nominal trajectory) evaluated at a given time
    % time
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to evaluate F_tilde           sec
    %
    % OUTPUT:       
    %    
    %  f_tilde_out      - taylor expanded dynamics matrix for S/C        (4x4)
    %                                     
    % Coupling:
    %   
    %  None 
    %
    %==========================================================================
    
    f_tilde_out = eye(4) + C.delta_t*[0, 1, 0, 0;
                                      C.mu * (3 * (cos(C.n * t))^2 - 1) / (C.r0^3), 0, ...
                                      3 * C.mu * sin(2 * C.n * t)/(2 * C.r0^3), 0;
                                      0, 0, 0, 1;
                                      3 * C.mu * sin(2 * C.n * t) / (2 * C.r0^3), 0, ...
                                      C.mu * ( 3 * (sin(C.n * t))^2 - 1) / (C.r0^3), 0];
end
  
%% Observation Matrix 
 function [h_matrix] = Build_Full_H_Matrix(C, time)
    %==========================================================================
    % [h_matrix] = Build_Full_H_Matrix(C, time)
    %
    % Consructs the observation matrix for all of the stations at a given
    % time
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %
    % OUTPUT:       
    %    
    %  h_matrix      - matrix of gps satellite orbit parameters         (nx25)
    %                                     
    % Coupling:
    %   
    % H_i => Uses this to construct the observation for each individual ground 
    %        station  
    %
    %==========================================================================
    
    h_matrix = [];
    for i = 1:12
        if Check_Valid_Pass(C, time, i)
            h_matrix = [h_matrix; H_i(C, time, i)];
        else
            h_matrix = [h_matrix; NaN(3, 4)];
        end
    end   
 end 
 
 
 function bool = Check_Valid_Pass(C, t, i)
    %==========================================================================
    % [bool] = check_valid_pass(C, t, i)
    %
    % Determines whether a spacecraft is visible within a given stations
    % viewing cone
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  t                   - time at which to initialize h matrix        sec
    %  i                   - number of station to generate h matrix for  int [1-12] 
    %
    % OUTPUT:       
    %    
    %  bool               - boolean as to whether S/C is visible to stn  (boolean)
    %                                     
    % Coupling:
    %
    %  phi_i              - finds elevation angle of spacecraft as viewed
    %                       from station
    %  theta_i            - finds angle of station in inertial frame
    %
    %==========================================================================
    
    phi = Phi_I(C, t, i); 
    theta = Theta_I(C, t, i); 
    min = -pi / 2 + theta;
    max = pi / 2 + theta;
    
    if min <= phi <= max
        bool = true;
    else 
        bool = false;
    end 
 end

 
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
    
    H_out = [(C.r0 * cos(C.n * t) - X_i(C, t, i)) / Rho_Nom(C, t, i), 0, ... 
            (C.r0 * sin(C.n * t) - Y_i(C, t, i)) / Rho_Nom(C, t, i), 0; 
           -(C.r0 * sin(C.n * t) - Y_i(C, t, i)) * Q_Nom(C, t, i) / Rho_Nom(C, t, i)^3, ...
            (C.r0 * sin(C.n * t) - X_i(C, t, i))/ Rho_Nom(C, t, i), ...
            (C.r0 * cos(C.n * t) - X_i(C, t, i)) * Q_Nom(C, t, i) / Rho_Nom(C, t, i)^3, ...
            (C.r0 * sin(C.n * t) - Y_i(C, t, i)) / Rho_Nom(C, t, i);
           -(C.r0 * sin(C.n * t) - Y_i(C, t, i)) / Rho_Nom(C, t, i)^2, 0, ...
            (C.r0 * cos(C.n * t) - X_i(C, t, i)) / Rho_Nom(C, t, i)^2, 0];
end    


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
    rho = sqrt(C.Re^2 + C.r0^2 - 2 * C.r0 * (X_i(C,t,i) * cos(C.n *t) + ...
               Y_i(C, t, i) * sin(C.n * t)));
end


function q = Q_Nom(C, t, i) 
    %==========================================================================
    % [q] = Q_nom(C, time)
    %
    % Solves a knarly bit of math necessary to generate the H matrix
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
    
    q = (C.r0 * cos(C.n * t) - X_i(C, t, i)) * (C.n * C.r0 * cos(C.n * t) - dY_i(C, t, i)) - ...
        (C.r0 * sin(C.n * t) - Y_i(C, t, i)) * (-C.n * C.r0 * sin(C.n * t) - dX_i(C, t, i));
end

 
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
    
    phi = atan2((C.r0 * sin(C.n * t) - Y_i(C, t, i)) , (C.r0 * cos(C.n * t) - X_i(C, t, i)));
 end
 
 
 %% Grounds Station Positions and angles
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








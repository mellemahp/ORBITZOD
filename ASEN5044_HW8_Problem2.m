
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
%C.Q_true = 
                                                   
%% Generate Truth States

% create the time vector for two periods of the orbit
times = 0:C.delta_t:14000; %C.period * 2;

% initialize the initial conditions of the orbit and perturbation vector
state = [C.r0; 0; 0; C.r0 * sqrt(C.mu / C.r0^3)];
p_vec= [0., 0.075, 0, -0.021]';

% set up nominal trajectory 
states_nom = [C.r0 * cos(C.n * times); 
              -C.r0 * C.n * sin(C.n * times);
              C.r0 * sin(C.n * times); 
              C.r0 * C.n * sin(C.n * times)];

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

% Construct perturbation measurements
msrs_p= [];
for t_index = 1:(length(times)) 
    h_mat_p = Build_Full_H_Matrix_Delta(C, times(t_index));
    msrs_p(:, t_index) = h_mat_p * p_vec(:, t_index);
end

% Construct nominal measurements 
msrs_nom = [];
for t_index = 1:(length(times))
    msrs_nom(:, t_index) = Build_Full_H_Matrix_Nom(C, times(t_index));
end


%% Plot Measurements
figure()

x = 1;
for i = 1:12
    subplot(3, 1, 1)
    hold on 
    scatter(times, msrs_nom(x,:) + msrs_p(x,:))
    ylabel('Range (km)')
    subplot(3, 1, 2) 
    hold on 
    scatter(times, msrs_nom(x + 1,:) + msrs_p(x + 1,:))
    ylabel('Range Rate (km/s)')
    subplot(3, 1, 3)
    hold on
    scatter(times, msrs_nom(x + 2,:) + msrs_p(x + 2,:))
    ylabel('Phi (rad)')
    xlabel('Time (sec)')
    x = x + 3; 
end

x = 2;
                    

%% Truth Modeling (Dynamics) 
 
load('orbitdeterm_finalproj_KFdata.mat')
x_true = [6678, 0, 0, C.r0 * sqrt(C.mu / C.r0^3)];
x_true = x_true + p_vec(:,1)';
t_span = [0, 10];
for k = 1:length(times)-1
    w_k = mvnrnd([0, 0], Qtrue, 1);
    [out_times, out_states] = ode45(@(t, x) Full_Nonlinear_Dynamics(C, t, x, w_k), t_span, x_true(k, :));
    x_true(k + 1, :) = out_states(end, :);
end

%% Plot dynamics truth model

figure()
y_labels = ["X (km)", "Xdot (km/s)", "Y (km)", "Ydot (km/s)"];
for i = 1:4
    subplot(4,1,i)
    plot(times, x_true(:,i))
    ylabel(y_labels(i))
    xlabel("Time (sec)")
end

%% Observation Truth Model 

msrs_true = []; 
for k = 1:length(times)-1 
    msrs_true(:, k+1) = Get_Msrs_True(C, x_true(k+1,:), times(k+1), Rtrue);
end

%% Plot Truth Measurements

figure()

x = 1;
for i = 1:12
    subplot(3, 1, 1)
    hold on 
    scatter(times, msrs_true(x,1:length(times)))
    ylabel('Range (km)')
    subplot(3, 1, 2) 
    hold on 
    scatter(times, msrs_true(x + 1,1:length(times)))
    ylabel('Range Rate (km/s)')
    subplot(3, 1, 3)
    hold on
    scatter(times, msrs_true(x + 2, 1:length(times)))
    ylabel('Phi (rad)')
    xlabel('Time (sec)')
    x = x + 3; 
end

%% FUNCTIONS
% ====================================================================================
P_0 = 1000 * eye(4); 
Q = eye(4);

%% pull out msrs

%% Run Kalman Filter 
msrs_true = msrs_true(:,2:end);
msrs_true_useful_vec = [];
H_block = [];
R_block=[];
xp = state;
P=P_0;
for t = 1:length(msrs_true)
    msrs_true_useful_vec = [];
    H_block = []; 
    R_block = [];
    for stn_num = 1:12
        j = 1:3:36;
        i = j(stn_num);
        if ~isnan(msrs_true(i,t))
            msrs_true_useful_vec = [msrs_true_useful_vec;msrs_true(i:i+2,t)];
            H_block = [H_block; H_i(C, times(t), stn_num)];
            R_block = blkdiag(R_block,Rtrue); 
        end
    end
    if ~isempty(msrs_true_useful_vec)
        F= F_tilde(C,t);
        [P(:,:,t+1), xp(:,t+1)] = Kalman_Step(F, H_block, xp(:,t), P(:,:,t), msrs_true_useful_vec, Q, R_block);
    end
end 


%% Kalman Filter Test 



[x_p, P, K] = Kalman_Filter2(C, state, P_0, Q, msrs_true, times, msrs_nom);

%% Kalman Filter 

function [xp, P, k] = Kalman_Filter2(C, xp_0, P_0, Q, data, times, msrs_nom)  
    % Preinitialization of cov, mean
    xp = xp_0;
    P = P_0;
    K = zeros(size(xp_0,1), size(msrs_nom, 1));

    for t = 1:(length(times) - 1)
        F_t = F_tilde(C, times(t));
        H_t = Build_Full_H_Matrix_Delta(C, times(t));
        RESID = msrs_nom(:, t) - data(:, t); 
        [P(:, :, t + 1), K(:, :, t + 1), xp(:, t + 1)] = Kalman_Step(F_t, H_t, xp(:, t), P(:, :, t), RESID, Q); 
    end
end


function [P, xp] = Kalman_Step(F, H, xp, P, msr, Q,R)
    % Prediction Step 
    xm = F * xp;
    Pm = F * P * F' + Q;
    K = Pm * H' * inv(H * Pm * H' + R);
    
    % Correction Step 
    xp = xm + K * (msr - H * xm);
    P = (eye(4) - K * H) * Pm;   
end


%% Truth Obs 

function [h_matrix] = Get_Msrs_True(C, sc_state, time, R) 
    h_matrix = [];
    for i = 1:12
        if Check_Valid_Pass_Truth(C, sc_state, time, i)
            h_new = True_Msrs(C, sc_state, time, i) + mvnrnd([0, 0, 0], R, 1)';
            h_matrix = [h_matrix; h_new];
        else
            h_matrix = [h_matrix; NaN(3, 1)];
        end
    end   
end 


 function bool = Check_Valid_Pass_Truth(C, sc_state, t, i)
    u = [sc_state(1) - X_i(C, t, i), sc_state(3) - Y_i(C, t, i), 0];
    v = [X_i(C, t, i), Y_i(C, t, i), 0];
    
    ThetaBetween = atan2(norm(cross(u,v)),dot(u,v));
    
    if -pi/2 < ThetaBetween && ThetaBetween < pi/2
        bool = true;
    else 
        bool = false;
    end 
 end


function [msrs_out] = True_Msrs(C, sc_state, t, i)
    msrs_out = [sqrt((sc_state(1) - X_i(C, t, i))^2 + (sc_state(3) - Y_i(C, t, i))^2); 
                ((sc_state(1) - X_i(C, t, i)) * (sc_state(2) - dX_i(C, t, i)) + ... 
                (sc_state(3) - Y_i(C, t, i)) * (sc_state(4) - dY_i(C, t, i))) / ... 
                sqrt((sc_state(1) - X_i(C, t, i))^2 + (sc_state(3) - Y_i(C, t, i))^2);
                atan2(sc_state(3) - Y_i(C, t, i), sc_state(1) - X_i(C, t, i))];
end


%% Full Nonlinear Dynamics 

function [state_acc] = Full_Nonlinear_Dynamics(C, t, state_m, noise) 
    % Generates accelarations based on full nonlinear EOM
    state_acc(1) = state_m(2);
    state_acc(2) = -C.mu * state_m(1) / sqrt(state_m(1)^2 + state_m(3)^2)^3 + noise(1);
    state_acc(3) = state_m(4);
    state_acc(4) = -C.mu * state_m(3) / sqrt(state_m(1)^2 + state_m(3)^2)^3 + noise(2);
    state_acc = state_acc';
end

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
    % OUTPUT:       zeros(size(xp_0,1), size(msrs_nom,1))'; 
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

%% Obs Matrix Nominal 

function [h_matrix] = Build_Full_H_Matrix_Nom(C, time) 
    %==========================================================================
    % [h_matrix] = Build_Full_H_Matrix_Nom(C, time)
    %
    % Consructs the Nominal observation matrix for all of the stations at a given
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
            h_matrix = [h_matrix; Nom_Msrs(C, time, i)];
        else
            h_matrix = [h_matrix; NaN(3, 1)];
        end
    end   
end 


function [msrs_out] = Nom_Msrs(C, t, i)
    %==========================================================================
    % [msrs_out] = Nom_msrs(C, t, i)
    %
    %  Finds the nominal measurements at a given time for a given station
    % 
    % INPUT:               Description                                   Units
    %
    %  C                   - Constants Object                            NA
    %  time                - time at which to initialize h matrix        sec
    %  i                  - number of station to generate h matrix for   int [1-12] 
    %
    % OUTPUT:       
    %    
    %                  - Range from station to nominal trajectory     km
    %                                     
    % Coupling:
    %
    %  X_i                - Gives the X location of station i at a given
    %                       time
    %  Y_i                - Gives the Y location of station i at a given
    %                       time
    %
    %==========================================================================
    msrs_out = [Rho_Nom(C, t, i); ((X_Nom(C, t) - X_i(C, t, i)) * (X_Nom_Dot(C, t) - dX_i(C, t, i)) + ...
                (Y_Nom(C, t) - Y_i(C, t, i)) * (Y_Nom_Dot(C, t) - dY_i(C, t, i))) / Rho_Nom(C, t, i); 
                atan2((Y_Nom(C, t) - Y_i(C, t, i)),(X_Nom(C, t) - X_i(C, t, i)))];

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
    rho = sqrt((X_Nom(C, t) - X_i(C, t, i))^2 + (Y_Nom(C, t) - Y_i(C, t, i))^2);
end


function [x_nom_out] = X_Nom(C, t)
    % Calculates the nominal x position at a given time
    x_nom_out = C.r0 * cos(C.n * t); 
end


function [x_nom_dot_out] = X_Nom_Dot(C, t)
    % Calculates the nominal x velocity at a given time
    x_nom_dot_out = -C.r0 * C.n * sin(C.n * t); 
end


function [y_nom_out] = Y_Nom(C, t) 
    % Calculates the nominal y position at a given time
    y_nom_out = C.r0 * sin(C.n * t);
end 


function [y_nom_dot_out] = Y_Nom_Dot(C, t)
    % Calculates the nominal y velocity at a given time
    y_nom_dot_out = C.r0 * C.n * cos(C.n * t);
end 


%% Observation Matrix (Perturbation)
 function [h_matrix] = Build_Full_H_Matrix_Delta(C, time)
    %==========================================================================
    % [h_matrix] = Build_Full_H_Matrix(C, time)
    %
    % Consructs the perturbation observation matrix for all of the stations at a 
    % given time
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
    u = [X_Nom(C, t) - X_i(C, t, i), Y_Nom(C,t) - Y_i(C, t, i), 0];
    v = [X_i(C, t, i), Y_i(C, t, i), 0];
    
    ThetaBetween = atan2(norm(cross(u,v)),dot(u,v));
    
    if -pi/2 < ThetaBetween && ThetaBetween < pi/2
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
    H_out = [(X_Nom(C, t) - X_i(C, t, i)) / Rho_Nom(C, t, i), 0, ... 
             (Y_Nom(C, t) - Y_i(C, t, i)) / Rho_Nom(C, t, i), 0;   
             (Y_i(C, t, i) - Y_Nom(C, t)) * Q_Nom(C, t, i) / Rho_Nom(C, t, i)^3, ...  
             (X_Nom(C, t) - X_i(C, t, i))/ Rho_Nom(C, t, i), ...
             (X_Nom(C, t) - X_i(C, t, i)) * Q_Nom(C, t, i) / Rho_Nom(C, t, i)^3, ...  
             (Y_Nom(C, t) - Y_i(C, t, i)) / Rho_Nom(C, t, i);        
            -(Y_Nom(C, t) - Y_i(C, t, i)) / Rho_Nom(C, t, i)^2, 0, ...
             (X_Nom(C, t) - X_i(C, t, i)) / Rho_Nom(C, t, i)^2, 0];
end    


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








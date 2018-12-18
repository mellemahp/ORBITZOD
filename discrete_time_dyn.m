clc 
clear all
close all 

%% Constants
C = Constants(); 
                                                   
%% Generate Nominal States

% create the time vector for two periods of the orbit
times = 0:C.delta_t:14000;

% initialize the initial conditions of the orbit and perturbation vector
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

Plot_States(times, p_vec', "Perturbation")


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

Plot_Measurements(times, msrs_nom + msrs_p)
                    

%% Truth Modeling (Dynamics with NO NOISE) 
 
x_true = [6678, 0, 0, C.r0 * sqrt(C.mu / C.r0^3)];
x_true = x_true + p_vec(:,1)'; 
t_span = [0, 10];
for k = 1:length(times)-1
    w_k = [0, 0];
    [out_times, out_states] = ode45(@(t, x) Full_Nonlinear_Dynamics(C, t, x, w_k), t_span, x_true(k, :));
    x_true(k + 1, :) = out_states(end, :);
end

%% Plot dynamics truth model

Plot_States(times, x_true, "Truth")

%% Observation Truth Model (NO NOISE)

msrs_true = []; 
for k = 1:length(times)-1 
    msrs_true(:, k+1) = Get_Msrs_True(C, x_true(k+1,:), times(k+1), zeros(3,3));
end

%% Plot Truth Measurements

Plot_Measurements(times, msrs_true)



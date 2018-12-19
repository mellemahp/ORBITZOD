
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

%% Observation Truth Model 

msrs_true = []; 
for k = 1:length(times)-1 
    msrs_true(:, k+1) = Get_Msrs_True(C, x_true(k+1,:), times(k+1), Rtrue);
end

%% Run Kalman Filter

msrs_true = msrs_true(:,2:end);
P_0 = 100 * eye(4); 
xp = p_vec;
P= P_0;
Q = eye(4) * 10^-9;

[P, xp] = Kalman_Filter(C, times, msrs_true, xp, P, Rtrue, Q);

%% Plot Kalman Filter Results 

%figure()
%plot(states_nom(1, :) + xp(1, :), states_nom(3, :) + xp(3, :))
%hold on 
%plot(states_nom(1, :), states_nom(3, :))

%% Extended Kalman Filter

P_0 = diag(1e-3 * [10, 0.001, 10, 0.001]);  
Q = (eye(4) * 1.0e-9) / 1.1; 
istate = [C.r0, 0, 0, C.r0 * sqrt(C.mu / C.r0^3)]';

[xp, P, ey, S] = EKF(C, istate, P_0, times, msrs_true, Q, Rtrue);


%% NEES and NIS Test
ex = x_true' - xp
[ex, ey, mu_ex, mu_ey] = NEESnNIS(ex, ey, P, S, msrs_true)


%% Plot Extended Kalman Filter Results 

figure()
plot(xp(1, :), xp(3,:))
hold on 
plot(states_nom(1, :), states_nom(3, :))


%% extract data from ydata cell array 

msrs_corrected = Make_Data_Useful(ydata);

%% Use a linear Kalman Filter on Given Data

P_0 = diag(1e-3 * [10, 0.001, 10, 0.001]);  
xp = p_vec;
P= P_0;
Q = (eye(4) * 1.0e-9) / 1.1; 
msrs_corrected = msrs_corrected(:,2:end);

[P, xp] = Kalman_Filter(C, times, msrs_corrected, xp, P, Rtrue, Q);

figure()
plot(states_nom(1, :) + xp(1, :), states_nom(3, :) + xp(3, :))
hold on 
plot(states_nom(1, :), states_nom(3, :))


Plot_States(times, (states_nom + xp)', "Perturbation")




















 


















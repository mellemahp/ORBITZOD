% Unscented Kalman Filter 

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
% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;
%%



istate = [C.r0, 0, 0, C.r0 * sqrt(C.mu / C.r0^3)]';
% Step 2: Define noise assumptions
G = [0 0;1 0;0 0;0 1]; 
Omega = C.delta_t * G; 
Q = Qtrue / 1.1; 
Q_Om = Omega * Q * Omega'; 
R = diag([10000 0.0100]);

N=length(x_true)
T=times
% Step 3: Initialize state and covariance
x = zeros(4,N); % Initialize size of state estimate for all k
xt(:,1) = istate ; % Set initial state estimate
P0 =diag(1e-3 * [10, 0.001, 10, 0.001]); % Set initial error covariance

yt = zeros(2, N); % Initialize size of output vector for all k
for k = 2:N
F = F_tilde(C, times(k));
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
    [~, out_states] = ode45(@(t, x) Full_Nonlinear_Dynamics(C, t, x, w_k), t_span, xp(:, k)', opts);
end

%%
% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 2:N
% Step 1: Generate the sigma-points
sP = chol(P,'lower'); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
         x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+chi_p(2,i)*T;... //x1
             chi_p(2,i)+((chi_p(1,i)*(chi_p(4,i))^2)-((G*M)/(chi_p(1,i))^2))*T;...//x2
             chi_p(3,i)+chi_p(4,i)*T;... //x3
             chi_p(4,i)-((2*chi_p(4,i)*chi_p(2,i))/chi_p(1,i))*T];   %//x4
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:)); ...
         (chi_m(3,:))];
y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end


    
%%
 %Plots
figure()
plot(xt(1, :), xt(3,:))
hold on 
plot(states_nom(1, :), states_nom(3, :))


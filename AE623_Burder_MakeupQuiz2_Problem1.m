%%AE623 - Makeup Quiz Problem 1
%Divinaa Burder

%% Before filter execution
close all
clc
clear all

% System properties
T = 0.1; % Sampling time
N = 100; % Number of time steps for filter

% Step 2: Define noise assumptions
Q_true = 1e-10;
Q_est = 1e-10;

% Step 3: Initialize state and covariance
xhat = zeros(2, N); % Initialize size of state estimate for all k
xhat(:,1) = [0;0]; % Set initial state estimate
P0 = eye(2); % Set initial error covariance

% Also calculate measurement vector
w = sqrt(Q_true)*randn(1, N); % Generate random process noise (from assumed Q)

xt = zeros(1, N); % Initialize size of true state for all k
xt(:,1) = 1; % Set true initial state
yt = zeros(1, N); % Initialize size of output vector for all k
for k = 2:N
xt(:,k) = 0.9*xt(:,k-1) + w(:,k-1);
yt(:,k) = xt(:,k-1);
end

%% Initialize and run EKF for comparison
P = P0;
 % Linear prediction
for k = 2:N
A = [xhat(2,k-1) xhat(1,k-1);0 1];
% Prediction
x_m = A*xhat(:,k-1)+ w(:,k-1);
P_m = A*P*A' + Q_est;
% Observation
y_m = xhat(1,k-1);
C = [1 1];
% Measurement Update
K = P_m*C'/(C*P_m*C'); % Calculate Kalman gain
xhat(:,k) = xhat(:,k-1)+K*(yt(:,k)-y_m);
P = (eye(2)-K*C)*P_m; % Update covariance estimate
end


%% Display results
figure(1);
t = T*(1:N);

figure(1)
plot(t,xhat(2,:),'b-', t,0.9*ones(1,length(xhat)),'r--', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('\phi (rad)')
grid on
legend('EKF','True')
title('AE623 - Burder - Makeup Quiz Problem 1 - \phi _E_K_F vs \phi _t_r_u_e')

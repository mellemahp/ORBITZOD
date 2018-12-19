
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

figure()
plot(states_nom(1, :) + xp(1, :), states_nom(3, :) + xp(3, :))
hold on 
plot(states_nom(1, :), states_nom(3, :))

%% Extended Kalman Filter

P_0 = diag(1e-3 * [10, 0.001, 10, 0.001]);
G = [0 0;1 0;0 0;0 1]; 
Omega = C.delta_t * G; 
Q = Qtrue / 1.1; 
Q_Om = Omega * Q * Omega'; 
istate = [C.r0, 0, 0, C.r0 * sqrt(C.mu / C.r0^3)]';

[xp, P, ey, S] = EKF(C, istate, P_0, times, msrs_true, Q_Om, Rtrue);



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

%% NEES and NIS Test for Linear Kalman Test

clc 
clear all
close all

%[epsilonx, epsilony, mu_ex, mu_ey] = NEESnNIS(ex, ey, P, S, msrs_true)

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


%% Run Kalman Filter


%msrs_true = msrs_true(:,2:end);

for j = 1:30
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
P_0 = blkdiag([10000],[100],[10000],[100]); 
xp = p_vec;
P= P_0;
G = [0 0;1 0;0 0;0 1]; 
Omega = C.delta_t * G; 
Q = Qtrue / 1.1; 
Q_Om = Omega * Q * Omega'; 
[P, xp, S, H, msrs] = Kalman_Filter(C, times, msrs_true(:,2:end), xp, P, Rtrue, Q_Om);

ex = x_true' - (xp + states_nom);
yHat = (H * xp);
ey = msrs(1:3,:) - yHat(:,2:end);

    for i = 1:201
        % NEES Test
        epsilon_x(:,i) = ex(:,i)' * inv(P(:,:,i)) * ex(:,i);
        mu_epsilonx(:,i) = mean(epsilon_x(:,i));
        
        %NIS Test 
        epsilon_y(:,i) = ey(:,i)' * inv(S{i}) * ey(:,i);
        mu_epsilony(:,i) = mean(epsilon_y(:,i));
    end
end

r1 = chi2inv(0.05/2,4*j)/j
r2 = chi2inv(1-0.05/2,4*j)/j

figure()
scatter(times(1:201),mu_epsilonx)
hold on
plot(times(1:201),r1*ones(1,201),'--',times(1:201),r2*ones(1,201),'--')
title('NEES Test')

figure()
scatter(times(1:201),mu_epsilony)
hold on
plot(times(1:201),r1*ones(1,201),'--',times(1:201),r2*ones(1,201),'--')
title('NIS Test')
    



%% NEES and NIS Test for Extended Kalman Filter
clear all
close all
clc

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

load('orbitdeterm_finalproj_KFdata.mat')
p_vec= [0., 0.075, 0, -0.021];
x_true = [6678, 0, 0, C.r0 * sqrt(C.mu / C.r0^3)];
x_true = x_true + p_vec;
t_span = [0, 10];
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
for j = 1:10
    
    for k = 1:length(times)-1
        w_k = mvnrnd([0, 0], Qtrue, 1);
        [out_times, out_states] = ode45(@(t, x) Full_Nonlinear_Dynamics(C, t, x, w_k), t_span, x_true(k, :), opts);
        x_true(k + 1, :) = out_states(end, :);
    end
    msrs_true = []; 

    for k = 1:length(times)-1 
        msrs_true(:, k+1) = Get_Msrs_True(C, x_true(k+1,:), times(k+1), Rtrue);
    end

    % Extended Kalman Filter

    P_0 = diag(1e-3 * [10, 0.001, 10, 0.001]);
    G = [0 0;1 0;0 0;0 1]; 
    Omega = C.delta_t * G; 
    Q = Qtrue; 
    Q_Om = Omega * Q * Omega'; 
    istate = x_true(1,:)';
    P = P_0;

    [xp_mc, P, ey, S] = EKF(C, istate, P_0, times, msrs_true, Q_Om, Rtrue);
    
    ex = x_true' - xp_mc;
    ex_log(:,:,j) = ex
    ey_log(:,:,j) = ey
end


%% NEES and NIS Test

[epsilonx, epsilony, mu_ex, mu_ey] = NEESnNIS(ex, ey, P, S, msrs_true)













 


















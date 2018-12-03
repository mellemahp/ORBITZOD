
clc
clear all
close all

%% Constants Definition
delta_t = 0.5;
time = linspace(0,300,300);
Omega_a = 0.045;
Omega_b = -0.045;

qw = 10;
W = qw*[2 0.05;0.05 0.5];

GamA = [0, 0; 1, 0; 0, 0; 0, 1];
GamB = [0, 0; 1, 0; 0, 0; 0, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_mat = @(omega) [0 1 0 0;0 0 0 -omega;0 0 0 1;0 omega 0 0];
A_a = A_mat(Omega_a)
A_b = A_mat(Omega_b)

F_mat = @(omega, dt) [1, sin(omega * dt) / omega, 0, -(1 - cos(omega * dt)) / omega;
                      0, cos(omega * dt), 0, -sin(omega * dt);
                      0, (1 - cos(omega * dt)) / omega, 1, sin(omega * dt) / omega;
                      0, sin(omega * dt), 0, cos(omega * dt)];


Fa = F_mat(Omega_a, delta_t)
Fb = F_mat(Omega_b, delta_t)


Z = @(dt, A, gamma, W_mat) dt * [-A, gamma * W_mat * gamma'; zeros(4, 4), A'];

Za = Z(delta_t, A_a, GamA, W);
Zb = Z(delta_t, A_b, GamB, W);

e_za = expm(Za);
e_zb = expm(Zb);

FinvQA = e_za(1:4, 5:8);
Qa = (Fa')' * FinvQA

FinvQB = e_zb(1:4, 5:8);
Qb = (Fb')' * FinvQB


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART B - i
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set random number seed for consistency between runs
rng(100)

Ha=[1, 0, 0, 0; 0, 0, 1, 0];
Ra=[20, 0.05; 0.05, 20];

load('hw8problem1_data')
xa = xasingle_truth;

% generate noisy measurements 
for i = 1:length(xasingle_truth)
    noise = mvnrnd([0 0], Ra);
    y_k(:,i) = Ha * xasingle_truth(:,i) + noise';
end

figure(1)
time = linspace(0, 20, 20);
plot(time,y_k(:, 1:20))

xlabel('Time (s)')
ylabel('y_a(k)')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART B - ii 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial mean and covariance guesses 
mu_a = [0 85 * cos(pi / 4), 0, -85 * sin(pi / 4)]';
P0_a = 900 * diag([10, 2, 10, 2]);

[xpa, Pa, Ka] = Kalman_Filter(mu_a, P0_a, Fa, Qa, y_k, Ha, Ra);

%%

est_state_error = xasingle_truth - xpa;

sigma1 = reshape(sqrt(Pa(1,1,:)),1,201);
sigma2 = reshape(sqrt(Pa(2,2,:)),1,201);
sigma3 = reshape(sqrt(Pa(3,3,:)),1,201);
sigma4 = reshape(sqrt(Pa(4,4,:)),1,201);


upper1 = mu_a(1,:) + 2*sigma1;
lower1 = mu_a(1,:) - 2*sigma1;

upper2 = + 2*sigma2;
lower2 = - 2*sigma2;

upper3 = mu_a(3,:) + 2*sigma3;
lower3 = mu_a(3,:) - 2*sigma3;

upper4 = + 2*sigma4;
lower4 = - 2*sigma4;


figure()
plot(time,est_state_error(1,1:20),time,upper1(1:20),'--r',time,lower1(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - State 1')
grid on
title('Easting Pos')

figure()
plot(time,est_state_error(2,1:20),time,upper2(1:20),'--r',time,lower2(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - State 2')
grid on
title('Easting Vel')

figure()
plot(time,est_state_error(3,1:20),time,upper3(1:20),'--r',time,lower3(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - State 3')
grid on
title('Northing Pos')

figure()
plot(time,est_state_error(4,1:20),time,upper4(1:20),'--r',time,lower4(1:20),'--r')

xlabel('Time(s)')
ylabel('Estimation state error - State 4')
grid on
title('Northing Vel')


%%%%%%%%%%%%%%%%%%%
%% PART C -- Setup 
%%%%%%%%%%%%%%%%%%%
mu_a_0 = [0, 85 * cos(pi / 4), 0, -85 * sin(pi / 4)]; 
P_a_0 = 900 * diag([10, 2, 10, 2]); 
mu_b_0 = [3200, 85 * cos(pi / 4), 3200, -85 * sin(pi / 4)];
P_b_0 = 900 * diag([11, 4, 11, 4]); 

R_d = [10 , 0.15; 0.15, 10]; 

%%%%%%%%%%%%%%%%
%% PART C -- i
%%%%%%%%%%%%%%%%

for i = 1:length(xadouble_truth)
    noise = mvnrnd([0 0],Ra);
    noise2 = mvnrnd([0 0], R_d); 
    y_ap(:,i)=Ha*xadouble_truth(:,i)+noise';
    y_d(:,i) = Ha*xadouble_truth(:, i) - Ha*xbdouble_truth(:,i) + noise2';
end

y_sk = [y_ap;y_d];
x_sk = [xadouble_truth;xbdouble_truth];


F_block = blkdiag(Fa,Fb);
Q_block = blkdiag(Qa,Qb);
R_block = blkdiag(Ra,R_d);
H_block = [Ha zeros(size(Ha));Ha -Ha];
P_block = blkdiag(P_a_0,P_b_0);
mu_block = [mu_a_0, mu_b_0]';

[xpd, Pd, Kd] = Kalman_Filter(mu_block, P_block, F_block, Q_block, y_sk, H_block, R_block);

est_state_error_d = x_sk - xpd;

Pad = Pd(1:4,1:4,:);
Pbd = Pd(5:8,5:8,:);


sigma1_a = reshape(sqrt(Pad(1,1,:)),1,201);
sigma3_a = reshape(sqrt(Pad(3,3,:)),1,201);
sigma1_b = reshape(sqrt(Pbd(1,1,:)),1,201);
sigma3_b = reshape(sqrt(Pbd(3,3,:)),1,201);


upper1 =+ 2*sigma1_a;
lower1 =  - 2*sigma1_a;

upper2 =  + 2*sigma3_a;
lower2 =  - 2*sigma3_a;

upper3 =  + 2*sigma1_b;
lower3 =  - 2*sigma1_b;

upper4 =  + 2*sigma3_b;
lower4 = - 2*sigma3_b;


figure()
plot(time,est_state_error_d(1,1:20),time,upper1(1:20),'--r',time,lower1(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - Aircraft A - Easting')
grid on
title('Easting Pos')

figure()
plot(time,est_state_error_d(3,1:20),time,upper2(1:20),'--r',time,lower2(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - Aircraft A - Northing')
grid on
title('Northing Pos')

figure()
plot(time,est_state_error_d(5,1:20),time,upper3(1:20),'--r',time,lower3(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - Aircraft B - Easting')
grid on
title('Easting Pos')

figure()
plot(time,est_state_error_d(7,1:20),time,upper4(1:20),'--r',time,lower4(1:20),'--r')

xlabel('Time(s)')
ylabel('Estimation state error - Aircraft B - Northing')
grid on
title('Northing Pos')

%%%%%%%%%%%%%%%%
%% PART C -- ii 
%%%%%%%%%%%%%%%%

y_sk = y_d;
H_block = H_block(3:4,:);

[xpd, Pd, Kd] = Kalman_Filter(mu_block, P_block, F_block, Q_block, y_sk, H_block, R_d);

Pad = Pd(1:4,1:4,:);
Pbd = Pd(5:8,5:8,:);


sigma1_a = reshape(sqrt(Pad(1,1,:)),1,201);
sigma3_a = reshape(sqrt(Pad(3,3,:)),1,201);
sigma1_b = reshape(sqrt(Pbd(1,1,:)),1,201);
sigma3_b = reshape(sqrt(Pbd(3,3,:)),1,201);


upper1 =  + 2*sigma1_a;
lower1 =  - 2*sigma1_a;

upper2 =  + 2*sigma3_a;
lower2 =  - 2*sigma3_a;

upper3 =  + 2*sigma1_b;
lower3 =  - 2*sigma1_b;

upper4 =  + 2*sigma3_b;
lower4 =  - 2*sigma3_b;


figure()
plot(time,est_state_error_d(1,1:20),time,upper1(1:20),'--r',time,lower1(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - Aircraft A - Easting')
grid on
title('Easting Pos')

figure()
plot(time,est_state_error_d(3,1:20),time,upper2(1:20),'--r',time,lower2(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - Aircraft A - Northing')
grid on
title('Northing Pos')

figure()
plot(time,est_state_error_d(5,1:20),time,upper3(1:20),'--r',time,lower3(1:20),'--r')
xlabel('Time(s)')
ylabel('Estimation state error - Aircraft B - Easting')
grid on
title('Easting Pos')

figure()
plot(time,est_state_error_d(7,1:20),time,upper4(1:20),'--r',time,lower4(1:20),'--r')

xlabel('Time(s)')
ylabel('Estimation state error - Aircraft B - Northing')
grid on
title('Northing Pos')



%%%%%%%%%%%%%%%%%
%% PART C -- iii
%%%%%%%%%%%%%%%%%%





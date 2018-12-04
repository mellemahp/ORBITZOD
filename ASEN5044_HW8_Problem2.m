
clc 
clear all
close all 

%% Constants
C.mu = 3.986e5; % km^3/sec^2
C.r0 = 6678; % km 
C.period = 2 * pi * sqrt(C.r0^3 / C.mu); % sec
C.n = sqrt(C.mu / C.r0^3);
C.delta_t = 10; % sec 
C.Re = 6379; %km 
C.omega_e = 2 * pi / 86400; % earth rotation rate
C.station_angles = arrayfun(@(i) (i - 1) * pi / 6, 1:12);

time = 1:C.delta_t:C.period; 
                                                   
%% States
state = [C.r0; 0; 0; C.r0 * sqrt(C.mu / C.r0^3)];
x_nom = C.r0 * cos(C.n * time);
y_nom = C.r0 * sin(C.n * time); 
G_tilde=C.delta_t * eye(4,2);

p_vec= [0.01, 0.0001, 0.01, 0.0001]';

for i = 1:(length(time)-1)
    p_vec(:,i+1) = F_tilde(C, time(i)) * p_vec(:, i);
end 

%% Plots
plot(time, p_vec(2,:))
plot((p_vec(1,:) + x_nom), (p_vec(3,:) + y_nom))
axis equal

%% Run simulation!
msrs = [];
for i = 1:(length(time)) 
    h_mat = build_full_h_matrix(C, time(i));
    msrs(:, i) = h_mat * p_vec(:,i);
end

%% Plot Measurements

x = 1;
for i = 1:12
    subplot(12, 3, x)
    plot(time, msrs(x,:))
    title(strcat('Stn', num2str(i), 'Range'))
    subplot(12, 3, x + 1) 
    plot(time, msrs(x + 1,:))
    title(strcat('Stn', num2str(i), 'Range Rate'))
    subplot(12, 3, x + 2)
    plot(time, msrs(x + 1,:))
    title(strcat('Stn', num2str(i), 'Phi'))
    x = x + 3; 
end 


                    
%% Functions

function f_t = F_tilde(C, t)
    f_t = eye(4) + C.delta_t*[0, 1, 0, 0;
                              C.mu * (3 * (cos(C.n * t))^2 - 1) / (C.r0^3), 0, ...
                              3 * C.mu * sin(2 * C.n * t)/(2 * C.r0^3), 0;
                              0, 0, 0, 1;
                              3 * C.mu * sin(2 * C.n * t) / (2 * C.r0^3), 0, ...
                              C.mu * ( 3 * (sin(C.n * t))^2 - 1) / (C.r0^3), 0];
end
    

function x_i = X_i(C, t, i) 
    x_i = C.Re * cos(C.omega_e * t + C.station_angles(i));
end


function dx_i = dX_i(C, t, i) 
    dx_i = -C.Re * C.omega_e * sin(C.omega_e * t + C.station_angles(i));
end


function y_i = Y_i(C, t, i)
    y_i = C.Re * sin(C.omega_e * t + C.station_angles(i));
end


function dy_i = dY_i(C, t, i) 
    dy_i = C.Re * C.omega_e * cos(C.omega_e * t + C.station_angles(i));
end 


function theta = theta_i(C, t, i) 
    theta = atan(X_i(C, t, i) / Y_i(C, t, i));
end


function bool = check_valid_pass(C, t, i)
    phi = phi_i(C, t, i); 
    theta = theta_i(C, t, i); 
    min = -pi / 2 + theta;
    max = pi / 2 + theta;
    if min <= phi <= max
        bool = true;
    else 
        bool = false;
    end 
end


function rho = rho_nom(C, t, i) 
    rho = sqrt(C.Re^2 + C.r0^2 - 2 * C.r0 * (X_i(C,t,i) * cos(C.n *t) + ...
               Y_i(C, t, i) * sin(C.n * t)));
end


function q = Q_nom(C, t, i) 
    q = (C.r0 * cos(C.n * t) - X_i(C, t, i)) * (C.n * C.r0 * cos(C.n * t) - dY_i(C, t, i)) - ...
        (C.r0 * sin(C.n * t) - Y_i(C, t, i)) * (-C.n * C.r0 * sin(C.n * t) - dX_i(C, t, i));
end


function H = H_i(C, t, i) 
    H = [(C.r0 * cos(C.n * t) - X_i(C, t, i)) / rho_nom(C, t, i), 0, ... 
         (C.r0 * sin(C.n * t) - Y_i(C, t, i)) / rho_nom(C, t, i), 0; 
         -(C.r0 * sin(C.n * t) - Y_i(C, t, i)) * Q_nom(C, t, i) / rho_nom(C, t, i)^3, ...
         (C.r0 * sin(C.n * t) - X_i(C, t, i))/ rho_nom(C, t, i), ...
         (C.r0 * cos(C.n * t) - X_i(C, t, i)) * Q_nom(C, t, i) / rho_nom(C, t, i)^3, ...
         (C.r0 * sin(C.n * t) - Y_i(C, t, i)) / rho_nom(C, t, i);
         -(C.r0 * sin(C.n * t) - Y_i(C, t, i)) / rho_nom(C, t, i)^2, 0, ...
         (C.r0 * cos(C.n * t) - X_i(C, t, i)) / rho_nom(C, t, i)^2, 0];
end

 
 function phi = phi_i(C, t, i) 
    phi = atan2((C.r0 * sin(C.n * t) - Y_i(C, t, i)) , (C.r0 * cos(C.n * t) - X_i(C, t, i)));
 end
 
 function h = build_full_h_matrix(C, t)
    h = [];
    for i = 1:12
        if check_valid_pass(C, t, i)
            h = [h; H_i(C, t, i)];
        else
            h = [h; NaN(3, 4)];
        end
    end
 end 

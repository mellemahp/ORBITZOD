
clc 
clear all
close all 

%% Constants
mu=3.986e5;
r0 = 6678;
period = 2*pi*sqrt(r0^3/mu);
n = sqrt(mu/r0^3);
delta_t = 10;
time = 1:10:period;


I = eye(4);

F_tilde = @(t) eye(4) + delta_t*[0, 1, 0, 0;
                                 mu * (3 * (cos(n * t))^2 - 1) / (r0^3), 0, ...
                                 3 * mu * sin(2 * n * t)/(2 * r0^3), 0;
                                 0, 0, 0, 1;
                                 3 * mu * sin(2 * n * t) / (2 * r0^3), 0, ...
                                 mu * ( 3 * (sin(n * t))^2 - 1) / (r0^3), 0];

%% States
state = [r0; 0; 0; r0 * sqrt(mu / r0^3)];
x_nom = r0 * cos(n * time);
y_nom = r0 * sin(n * time); 
G_tilde=delta_t * eye(2);

p_vect = [.0, .075, .0, -0.021]';

for i = 1:length(time)
    state(:,i+1) = F_tilde(time(i)) * p_vect;   
end 


%% Plots
figure(1)
plot(x_nom + state(1,2:end),y_nom + state(3,2:end))
axis equal
figure()
plot(time, state(1,2:end))
figure()
plot(time, state(3,2:end)) 
size(time)







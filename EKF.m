function [X_out, P_out, eps_y_out, S_out] = EKF(Xm_p, Pm_p, Q, R, yHist, time, options)
%! Xm_p should be given as the initial state conditions (or guess at initial
%state condtiions
%! Pm_p should reflect confidence in Xm_p
%! Qm should be the noise matrix, initially. So far, this code doesn't
%reflect variant Q matrices.
%! R is process noise covariance
%! yHist is side-stacked vector of measurements
%! dt is constant time step.
%! um is constant or time history of forcing vector
numPoints = size(time,2);
numStates = 4;
numMeas = 3;
I = eye(numStates);
X_out = zeros(numStates, numPoints);
P_out = zeros(numStates, numStates, numPoints);
P_out(:,:,1) = Pm_p;
S_out = zeros(numMeas, numMeas, numPoints);
S_out(:,:,1) = R; %maybe a bad idea?
eps_y_out = zeros(numMeas, numPoints);

%extract measurements and stations
Lm = [0,0;10,0;0,0;0,10];
M = eye(numMeas);
%Qm = diag([10^(-3.5),10^(-6.8),10^(-3.5),10^(-6.8)])*10^(2.6);

for k = 1:numPoints
    
    y = yHist(1:3,k+1); %measurement for this time
    if k > 1
        dt = time(k) - time(k - 1);
        %Extract stuff from vectorized data first to keep code below simpler.
        %Do Kalman filter
        %Prediction
        Fm = make_F(Xm_p, dt);
        P_m = Fm * Pm_p * Fm' + Lm * Q * Lm';
        [~, X_m_hist] = ode45(@nonlinearOrbitODE, [0, dt], Xm_p, options); %i.e. our apriori state estimtate for x is a. Note that we propagate without noise.
        X_m = X_m_hist(end, :)';
    else
        dt = time(k);
        X_m = Xm_p;
        P_m = Pm_p;
    end

    %estimate corrections based on measurement
    stationObserving = yHist(4,k+1);
    if stationObserving ~= 0
        %function of th previous a posteriori and forcing last time.
        [H, yHat] = make_H_y(stationObserving, X_m, time(k)); %note that yHat is estimated without noise.
        S = (H * P_m * H' + M * R * M');
        K = P_m * H' / S;
        eps_y = y - yHat;
        Xm_p = X_m + K * eps_y; %i.e. innovation is measurement minus
        %predicted measurement.
        Pm_p = (I - K * H) * P_m;% * (I - K * H)' + K * R * K';
    else
        Xm_p = X_m;
        Pm_p = P_m;
    end
    

    %Write to output vectors
    X_out(:, k) = Xm_p;
    P_out(:,:,k) = Pm_p;
    eps_y_out(:,k) = eps_y;
    S_out(:,:,k) = S;
end


end


function STM = make_F(Xm_p, dt)

mu = 398600.;
x1 = Xm_p(1);
x3 = Xm_p(3);

STM = zeros(4, 4);

Atilde = zeros(4,4);
Atilde(1,:) = [0, 1, 0, 0];
Atilde(2,1) = -mu*(x1^2+x3^2)^(-3/2)+(3*mu*x1^2)*((x1^2+x3^2)^(-5/2));
Atilde(2,3) = (3*mu*x1*x3)*(x1^2+x3^2)^(-5/2);
Atilde(3,:) = [0, 0, 0, 1];
Atilde(4,1) = (3*mu*x1*x3)*(x1^2+x3^2)^(-5/2);
Atilde(4,3) = -mu*(x1^2+x3^2)^(-3/2) + (3*mu*x3^2)*((x1^2+x3^2)^(-5/2));

STM(:, :) = dt * Atilde + eye(4);
end

function [Htilde, y_est] = make_H_y(stationObserving, X, time)
%Constants
om_E = 2 * pi / 86400.;
R_E = 6378;

% Extract state variables from X

if stationObserving == 0
    Htilde = zeros(3,4);
    y_est = [0;0;0];
else
    theta_0 = (stationObserving - 1) * pi / 3.;
    x_s = R_E * cos(om_E * time + theta_0);
    y_s = R_E * sin(om_E * time + theta_0);
    x_dot_s = -R_E * om_E * sin(om_E * time + theta_0);
    y_dot_s = R_E * om_E * cos(om_E * time + theta_0);
    
    %Extract state variables
    x1 = X(1);
    x2 = X(2);
    x3 = X(3);
    x4 = X(4);
    
    f = (x1 - x_s) .* (x2 - x_dot_s) + (x3 - y_s) .* (x4 - y_dot_s);
    rho = sqrt((x1 - x_s).^2 + (x3 - y_s).^2);
    g = rho;
    
    df_dx1 = x2 - x_dot_s;
    df_dx2 = x1 - x_s;
    df_dx3 = x4 - y_dot_s;
    df_dx4 =  x3 - y_s;
    
    d_rho_dx_1 = (x1 - x_s) ./ sqrt((x1 - x_s).^2 + (x3 - y_s).^2);
    d_rho_dx_2 = 0.;
    d_rho_dx_3 = (x3 - y_s) ./ sqrt((x1 - x_s).^2 + (x3 - y_s).^2);
    d_rho_dx_4 = 0.;
    
    dphi_dx_1 = -(x3 - y_s) .* (x1 - x_s).^(-2) ./ (1 + ((x3 - y_s) ./ (x1 - x_s)).^2.);
    dphi_dx_3 = (1 ./ (x1 - x_s)) ./ (1 + ((x3 - y_s) ./ (x1 - x_s)).^2.);
    
    Htilde = zeros(3, 4);
    
    Htilde(1,1) = d_rho_dx_1;
    Htilde(1,3) = d_rho_dx_3;
    Htilde(2,1) = (df_dx1 * g- d_rho_dx_1*f) / (g^2);
    Htilde(2,2) = (df_dx2 * g - d_rho_dx_2*f) / (g^2);
    Htilde(2,3) = (df_dx3 * g - d_rho_dx_3*f) / (g^2);
    Htilde(2,4) = (df_dx4 * g - d_rho_dx_4*f) / (g^2);
    Htilde(3,1) = dphi_dx_1;
    Htilde(3,3) = dphi_dx_3;
    
    y_est = [rho;
        f / g;
        wrapToPi(atan2(x3 - y_s, x1 - x_s))]; %no noise because we are jus making our best guess.
end

end
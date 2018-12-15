function [xp, P, epyp, S_out] = EKF(mu_0, P_0,data, Q, R, L, tspan)

z = size(tspan,2);
xcount = length(mu_0);
ycount = 3;
I = eye(xcount);
xp = zeros(xcount, z);
P = zeros(xcount, xcount, z);
P(:,:,1) = P_0;
M = eye(ycount);
S_out = zeros(ycount, ycount, z);
S_out(:,:,1) = R; 
epyp = zeros(ycount, z);
mu = 398600;


for j = 1:z
    
    yk = data(1:3,j+1);
    if j > 1
        deltat = tspan(j) - tspan(j - 1);
        Fm = Ftilde(deltat);
        Pm = Fm * P_0 * Fm' + L * Q * L';
        tspan = [0 deltat]
        [t, statelog] = ode45(@NLOrbit,mu, tspan, mu_0); 
        xm = statelog(end, :)';
    else
        deltat = tspan(j);
        xm = mu_0;
        Pm = P_0;
    end

    Obs = data(4,j+1);
    if Obs ~= 0
        [yHat] = Build_Full_H_Matrix_Delta(Obs, tspan(j)); 
        S = (H * Pm * H' + M * R * M');
        K = Pm * H' / S;
        epy = yk - yHat;
        mu_0 = xm + K * epy;
        P_0 = (I - K * H) * Pm;
    else
        mu_0 = xm;
        P_0 = Pm;
    end
    
    xp(:, j) = mu_0;
    P(:,:,j) = P_0;
    epyp(:,j) = epy;
    S_out(:,:,j) = S;
end

end

function [xp, P, ey] = EKF(C, xp, P, times, msrs, Q, Rmat)
t_span = [0, 10];
z=size(times, 2);
ycount = 3;
ey = zeros(ycount, z);
w_k = [0, 0]; 

% Runs a Kalman Filter on a set of measurements
for k = 1:length(msrs)
    msrs_block = [];
    H_block = []; 
    R_block = [];
    for stn_num = 1:12
        j = 1:3:36;
        i = j(stn_num);
        if ~isnan(msrs(i, k))
            msrs_block = [msrs_block; msrs(i:i+2, k)]; 
            H_block = [H_block; H_i(C, times(k), stn_num)];
            R_block = blkdiag(R_block, Rmat); 
        end
    end
    F = F_tilde(C, times(k));
    % Prediction
    [~, out_states] = ode45(@(t, x) Full_Nonlinear_Dynamics(C, t, x, w_k), t_span, xp(:, k)'); 
    xm = out_states(end,:)';
    Pm = F * P(:, :, k) * F'+ Q;
    if ~isempty(msrs_block)
        % Correction
        S = (H_block * Pm * H_block' + R_block);
        K = Pm * H_block' / S;
        yHat = H_block * xm;
        e = msrs_block - yHat;
        xp(:, k+1) = xm + K * e;
        P(:,:,k+1) = (eye(4) - K * H_block) * Pm;
    else
        xp(:, k+1) = xm;
        P(:, :, k+1) = Pm;
    end
end

end

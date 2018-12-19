function [xp, P, ey] = EKF(C, xp, P, times, msrs, Q, Rmat)
t_span = [0, 10];
z = size(times, 2);
ycount = 3;
w_k = [0, 0]; 

% Runs a Kalman Filter on a set of measurements
for k = 1:length(msrs)
    pred_msrs = [];
    msrs_block = [];
    H_block = []; 
    R_block = [];
    
    % Prediction
    F = F_tilde(C, times(k));
    [~, out_states] = ode45(@(t, x) Full_Nonlinear_Dynamics(C, t, x, w_k), t_span, xp(:, k)); 
    xm = out_states(end,:)'; 
    Pm = F * P(:, :, k) * F'+ Q;
    
    for stn_num = 1:12
        j = 1:3:36;
        i = j(stn_num);
        if ~isnan(msrs(i, k))
            pred_msrs = [pred_msrs; True_Msrs(C, xm, times(k), stn_num)];
            msrs_block = [msrs_block; msrs(i:i+2, k)]; 
            H_block = [H_block; H_i(C, times(k), stn_num)];
            R_block = blkdiag(R_block, Rmat); 
        end
    end

    if ~isempty(msrs_block)
        % Correction
        S = (H_block * Pm * H_block' + R_block);
        K = Pm * H_block' / S;
        ey{:, k+1} = msrs_block - pred_msrs;
        xp(:, k+1) = xm + K * ey{:, k+1};
        P(:, :, k+1) = (eye(4) - K * H_block) * Pm;
    else
        xp(:, k+1) = xm;
        P(:, :, k+1) = Pm;
    end
end

end

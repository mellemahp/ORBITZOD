function [xp, P, K] = Kalman_Filter(mu_0, P_0, F, Q, data, H, R)
   
    % Preinitialization of cov, mean
    xp = mu_0;
    P = P_0;
    K = zeros(size(H * P_0))';

    for i = 1:(length(data) - 1)
        % Prediction Step 
        xm = F * xp(:, i);
        Pm = F * P(:,:,i) * F' + Q;
        K(:, :, i+1) = Pm * H' * inv(H * Pm * H' + R);
    
        % Correction Step 
        xp(:, i+1) = xm + K(:, :, i+1) * (data(:, i+1) - H * xm);
        P(:, :, i+1) = (eye(4) - K(:, :, i+1) * H) * Pm;
    end
end


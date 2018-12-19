function [P, xp] = Kalman_Step(F, H, xp, P, msr, Q,R)
    % Prediction Step 
    xm = F * xp;
    Pm = F * P * F' + Q;
    K = Pm * H' * inv(H * Pm * H' + R);
    
    % Correction Step 
    xp = xm + K * (msr - H * xm);
    P = (eye(4) - K * H) * Pm;   
end

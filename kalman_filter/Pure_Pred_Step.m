function [P, xp] = Pure_Pred_Step(F, xp, P, Q)
    % Runs pure prediction for an update
    xp = F * xp;
    P = F * P * F' + Q;
end
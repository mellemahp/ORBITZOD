function [P, xp, S_log, H_block, msrs_sansNaN] = Kalman_Filter(C, times, msrs, xp, P, Rmat, Q) 
    % Runs a Kalman Filter on a set of measurements
    for t = 1:length(msrs)
        residuals_bloc = [];
        H_block = []; 
        R_block = [];
        for stn_num = 1:12
            j = 1:3:36;
            i = j(stn_num);
            if ~isnan(msrs(i, t))
                resid_t = msrs(i:i + 2, t) - Nom_Msrs(C, times(t), stn_num);
                residuals_bloc = [residuals_bloc; resid_t]; 
                H_block = [H_block; H_i(C, times(t), stn_num)];
                R_block = blkdiag(R_block, Rmat); 
                msrs_sansNaN(:,t) = msrs(i:i + 2, t);
                
            end
        end
        F = F_tilde(C, t);
        if ~isempty(residuals_bloc)
            [P(:,:,t+1), xp(:,t+1)] = Kalman_Step(F, H_block, xp(:,t), P(:,:,t), residuals_bloc, Q, R_block);
            S = (H_block(1:3,1:4) * P(:,:,t) * H_block(1:3,1:4)' + R_block(1:3,1:3));
            S_log{t} = S;
        else 
            [P(:,:,t+1), xp(:,t+1)] = Pure_Pred_Step(F, xp(:,t), P(:,:,t), Q);
        end
    end
end


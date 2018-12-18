function [msrs_out] = True_Msrs(C, sc_state, t, i)
    msrs_out = [sqrt((sc_state(1) - X_i(C, t, i))^2 + (sc_state(3) - Y_i(C, t, i))^2); 
                ((sc_state(1) - X_i(C, t, i)) * (sc_state(2) - dX_i(C, t, i)) + ... 
                (sc_state(3) - Y_i(C, t, i)) * (sc_state(4) - dY_i(C, t, i))) / ... 
                sqrt((sc_state(1) - X_i(C, t, i))^2 + (sc_state(3) - Y_i(C, t, i))^2);
                atan2(sc_state(3) - Y_i(C, t, i), sc_state(1) - X_i(C, t, i))];
end


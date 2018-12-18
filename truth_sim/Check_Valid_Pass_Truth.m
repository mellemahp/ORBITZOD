function bool = Check_Valid_Pass_Truth(C, sc_state, t, i)
    u = [sc_state(1) - X_i(C, t, i), sc_state(3) - Y_i(C, t, i), 0];
    v = [X_i(C, t, i), Y_i(C, t, i), 0];
    
    ThetaBetween = atan2(norm(cross(u,v)),dot(u,v));
    
    if -pi/2 < ThetaBetween && ThetaBetween < pi/2
        bool = true;
    else 
        bool = false;
    end 
end


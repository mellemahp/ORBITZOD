function xdot = NLOrbit(mu,tspan, state)


xdot =  [state(2);- mu * state(1) / ((state(1)^2 + state(3)^2)^ (3/2));state(4);- mu * state(3) / ((state(1)^2 + state(3)^2)^ (3/2))];

end

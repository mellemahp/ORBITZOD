function Plot_States(times, states, prefix)
    % Plots the x, xdot, y, ydot states for a given dataset
    figure()
    
    titles = [" X (km)", " Xdot (km/s)", " Y (km)", " Ydot (km/s)"];
    for k = 1:4
        subplot(4,1,k)
        plot(times, states(:,k))
        title(strcat(prefix, titles(k)))
        xlabel("Time (sec)")
        if k == 2
            ylabel(strcat(prefix," Magnitude (km)               "))
        end 
    end
end


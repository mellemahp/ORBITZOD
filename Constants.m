classdef Constants
    % Initializes all necessary constants for the orbit determination
    % problem
    properties
        mu = 3.986e5; % km^3/sec^2
        r0 = 6678; % km 
        delta_t = 10; % sec 
        Re = 6378; %km 
        omega_e = 2 * pi / 86400; % earth rotation rate
        station_angles_0 = arrayfun(@(i) (i - 1) * pi / 6, 1:12);
    end
    methods 
        function p = period(self)
            p = 2 * pi * sqrt(self.r0^3 / self.mu); % sec
        end
        function nu = n(self)
            nu = 2 * pi / self.period;
        end
    end
end


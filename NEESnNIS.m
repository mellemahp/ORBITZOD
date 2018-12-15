function [epx, epy, mu_epx, mu_epy] = NEESnNIS(ex, ey, P, S, Obs)

count = size(P,4);
xcount = size(P,1);
ycount = size(S,1);
z = size(ex, 2);

epx = zeros(1, z, count);
mu_epx = zeros(1, z);
mu_epy = zeros(1, z);
for i = 1:count
   
    for s = 1:z
        epx(1, s, i) = (ex(:, s, i)' / P(:, :, s, i)) * ex(:, s, i); 
        epy(1, s, i) = ey(:, s, i)' / S(:, :, s, i) * ey(:, s, i);
    end
end

for s = 1:z
        mu_epx(1,s) = mean(epx(1,s,:));
        mu_epy(1,s) = mean(epy(1,s,:));
end

mu_epx(Obs == 0) = nan;

end

function [ex, ey, mu_ex, mu_ey] = NEESnNIS(ex, ey, P, S, Obs)

count = size(P,4);
xcount = size(P,1);
ycount = size(S,1);
z = size(ex, 2);

ex = zeros(1, z, count);
mu_ex = zeros(1, z);
mu_ey = zeros(1, z);
for i = 1:count
   
    for s = 1:z
        ex(1, s, i) = (ex(:, s, i)' / P(:, :, s, i)) * ex(:, s, i); 
        ey(1, s, i) = ey(:, s, i)' / S(:, :, s, i) * ey(:, s, i);
    end
end

for s = 1:z
        mu_ex(1,s) = mean(ex(1,s,:));
        mu_ey(1,s) = mean(ey(1,s,:));
end

mu_ex(Obs == 0) = nan;

end

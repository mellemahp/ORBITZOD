function [ex, ey, mu_ex, mu_ey] = NEESnNIS(ex, ey, P, S, msrs)

count = size(P,4);
xcount = size(P,1);
ycount = size(S,1);
z = size(ex, 2);

mu_ex = zeros(1, z);
mu_ey = zeros(1, z);
for i = 1:count
   
    epsilonx(:, i) = (ex(:, i)' * inv( P(:, :, i))) * ex(:,i); 
    epsilony(:, i) = ey{:, i}' * inv(S(:, :, i)) * ey{:, i};

end

for s = 1:z
        mu_ex(1,s) = mean(epsilonx(:,s));
        mu_ey(1,s) = mean(epsilony(:,s));
end

mu_ex(msrs == 0) = nan;

end

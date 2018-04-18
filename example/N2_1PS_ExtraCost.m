function [c, g] = N2_1PS_ExtraCost(p,idx,lambda)

p = reshape(p,[],3);
g = zeros(size(p));

idx = find(idx);

p = p(idx,:);
% s = std(p(:,2));
% mu = mean(p(:,2));

c = lambda*sum(p(:,1).*p(:,1));% + 100000*s;
g(idx,1) = 2*lambda*p(:,1);
% g(idx,2) = 100000/s*(p(:,2)-mu)/numel(idx);

g = g(:);

end
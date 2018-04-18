function [c, g] = H_ExtraCost(p,idx)

p = reshape(p,[],3);
g = zeros(size(p));

idx = find(idx);

p = p(idx,:);

% [~,idx2] = sort(p(:,2));
% p = p(idx2,:);
% idx = idx(idx2);

c = 5*sum(p(1:2,1).*p(1:2,1));% + 100000*s;
g(idx(1:2),1) = 10*p(1:2,1);
% g(idx,2) = 100000/s*(p(:,2)-mu)/numel(idx);

g = g(:);

end
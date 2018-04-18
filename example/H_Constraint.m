function [c,ceq,g,geq,initCorrections] = H_Constraint(p, idx)
% Ensure peaks next to each other are within some range of each other
% g : should be nx3 with the indices and values that can be passed to
% sparse -> sparse(g(:,1),g(:,2),g(:,3))
ceq = [];
geq = [];

% N = numel(A);
p = reshape(p,[],3);
A = p(:,1);
mu = p(:,2);

idx = find(idx);

A = A(idx);
mu = mu(idx);

% sort by wavelength
[~,idx3] = sort(mu);
A = A(idx3);
idx = idx(idx3); 

c(2) = A(3) - 20*A(2)*A(2);
c(1) = A(1) - 5*A(2)*A(2);
c = c(:);
i = [1,1,2,2];
j = idx([3,2,1,2]);
v = [1,-20*A(2),1,-5*A(2)];
g = [j(:),i(:),v(:)];

if nargout > 4
    initCorrections = zeros(3,1);    
    initCorrections(idx([1,3])) = max(c,0);
end

end
function [c,ceq,g,geq,initCorrections] = N2_2PS_Constraint(p, delta, idx)
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

% Find main peak at ~336;
[d,idx336] = min(abs(mu-336.7));
if d > 2
    % line is not present, do not apply any constraints
    c = [];
    g = [];
    return
end

% Extract N2 peaks over 365 and less than 330. These peaks should all be less than 0.5*N2_336.
idx2 = find(mu > 365);
idx3 = find(mu < 330);
index = [idx336;idx2;idx3];
A = A(index);
mu = mu(index);
idx2 = idx(index);

% [~,idx3] = sort(mu);
% A = A(idx3);
% idx = idx(idx3); 

n = numel(A);
A = A(:);
% constraints
% 1 :  A(2:end) < (1-delta)*A(1)

c = A(2:n) - (1-delta)*A(1);

% gradient (as sparse matrix)



i = (1:(n-1));
i = [i, i];

j = idx2([ones(1,n-1), 2:n]);

v = [-(1-delta)*ones(1,n-1), ones(1,n-1)];
% g = sparse(j,i,[v1,v2,v2,v1], N, max(i));

g = [j(:),i(:),v(:)];

if nargout > 4
    initCorrections = zeros(size(idx));
    initCorrections(index(2:end)) = c;
end

end
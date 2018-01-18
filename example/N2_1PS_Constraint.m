function [c,ceq,g,geq] = N2_1PS_Constraint(p, delta, idx)
% Ensure peaks next to each other are within some range of each other
% g : should be nx3 with the indices and values that can be passed to
% sparse -> sparse(g(:,1),g(:,2),g(:,3))
ceq = [];
geq = [];

% N = numel(A);
p = reshape(p,[],3);
A = p(:,1); % peak amplitudes
mu = p(:,2); % peak positions

idx = find(idx);

% amplitude and centers for lines of this constraint
A = A(idx);
mu = mu(idx);
n = numel(A);


% constraints - ensure neighboring peaks are within some fraction of each
% other.
% 1 :  A(i) > (1-delta)*A(i+1)
% 2 :  A(i) < (1+delta)*A(i+1)
% 3 :  A(i+1) > (1-delta)*A(i)
% 4 :  A(i+1) < (1+delta)*A(i)

c = [(1-delta)*A(2:n) - A(1:n-1); 
     A(1:n-1) - (1+delta)*A(2:n);
     (1-delta)*A(1:n-1) - A(2:n);
     A(2:n) - (1+delta)*A(1:n-1)];

% gradient (as sparse matrix)
ci = (1:4*(n-1));
xi = 1:n-1;
xi = [xi, xi, xi, xi];% + find(idx,1,'first')-1;


ci = [ci, ci];
xi = idx([xi, xi+1]);
v1 = [-1*ones(1,n-1), ones(1,n-1)];
v2 = [(1-delta)*ones(1,n-1), -(1+delta)*ones(1,n-1)];
v = [v1,v2,v2,v2];
% g = sparse(j,i,[v1,v2,v2,v1], N, max(i));

g = [xi(:),ci(:),v(:)];


% additional contraints - the N2 1PS system has 3 primary bands in my
% system. The first and second bands each have a peak that should be larger
% than all others in the band. Constrain the problem to make this sure.

% [d,idx605] = min(abs(mu - 605.7));
% if d < 2
%     % indices of this band
%     idx2 = [idx605; find((mu < 603 | mu > 607) & mu < 620)];
%     
%     if numel(idx2) > 1
%         a = A(idx2);
%         idx2 = idx(idx2);
%         n2 = numel(a);
% 
%         c2 = a(2:end) - a(1);
% 
%         ci = 1:(n2-1);
%         ci = [ci, ci] + size(c,1);
% 
%         xi = idx2([ones(1,n2-1), 2:n2]);
%         v = [-ones(1,n2-1), ones(1,n2-1)];
% 
%         c = [c;  c2];
%         g = [ g; xi(:), ci(:), v(:)];
%     end
% end
% 
% [d,idx661] = min(abs(mu - 661));
% if d < 2
%     % indices of this band
%     idx2 = [idx661; find((mu < 659 | mu > 663) & mu < 690 & mu > 620)]; 
%     
%     if numel(idx2) > 1
%         a = A(idx2);
%         idx2 = idx(idx2);
%         n2 = numel(a);
% 
%         c2 = a(2:end) - a(1);
% 
%         ci = 1:(n2-1);
%         ci = [ci, ci] + size(c,1);
% 
%         xi = idx2([ones(1,n2-1), 2:n2]);
%         v = [-ones(1,n2-1), ones(1,n2-1)];
% 
%         c = [c;  c2];
%         g = [ g; xi(:), ci(:), v(:)];
%     end
% end



delta = 0.0;

[d,ip] = min(abs(mu - 605.7));
if d < 2
    ind = find(mu < 620);
    ng = numel(ind);
    if ng>1
        v = [ones(ng-1,1), -(1+delta)*ones(ng-1,1)];
        v(ind(1:ng-1)>=ip,:) = v(ind(1:ng-1)>=ip,[2,1]);
        
        ci = 1:(ng-1);
        ci = [ci, ci];
        xi = [1:ng-1, 2:ng];
        
%         full(sparse(ci,xi,v(:)))
%         A(ind)        
        c2 = full(sparse(ci,xi,v(:))*A(ind));
        
        ci = ci + size(c,1);
        xi = idx(ind(xi));
        
        c = [c; c2];
        g = [g; xi(:), ci(:), v(:)];
    end
end

% error('some err')

[d,ip] = min(abs(mu - 661));
if d < 2
    ind = find(mu > 620 & mu < 690);
    ng = numel(ind);
    if ng>1
        v = [ones(ng-1,1), -(1+delta)*ones(ng-1,1)];
        v(ind(1:ng-1)>=ip,:) = v(ind(1:ng-1)>=ip,[2,1]);
        
        ci = 1:(ng-1);
        ci = [ci, ci];
        xi = [1:ng-1, 2:ng];
        
        c2 = full(sparse(ci,xi,v(:))*A(ind));
        
        ci = ci + size(c,1);
        xi = idx(ind(xi));
        
        c = [c; c2];
        g = [g; xi(:), ci(:), v(:)];
    end
end


end
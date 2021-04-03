function X = mypinv(A,tol)
%PINV   Pseudoinverse.


[U,S,V] = svdecon(A);
s = diag(S);
if nargin < 2 
    tol = max(size(A)) * eps(norm(s,inf));
end
r1 = sum(s > tol)+1;
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
s = 1./s(:);
X = bsxfun(@times,V,s.')*U';

function [x,W,iter] = activeSet(P,q,A,b,x0)

% Active Set Method for a Inequality-Constrained Quadratic Problem
%
% [x,W,iter] = activeSet(P,q,A,b,x0)
%
% Quadratic Problem: min 1/2(x+delta)'P(x+delta) + q'(x+delta) s.t. A(x+delta) <= b
%
% Input arguments: 
% P is n-by-n symmetric, positive definite matrix
% q is n-by-1 column vector
% A is m-by-n matrix
% b is m-by-1 column vector
% x0 is feasible starting point (n-by-1)
% where n is number of optimized variables and m is number of inequality constraints
%
% Output arguments:
% x is optimal solution
% W is final set of active constraints
% iter is number of iterations

%%

if size(P,1) ~= size(P,2)
  error('P must be square matrix')
end

eiv = eig(P);

if min(eiv) <= 0
  error('P must be positive definite matrix')
end

%SIZE CHECK

if size(P,1) ~= length(q)
  error('P and q must be matrices with the same number of rows.')
end

%Unconstrained Quadratic Programming

if isempty(A) == true && isempty(b) == true && isempty(x0) == true
  x = -inv(P)*q; W = []; iter = [];
  return
end

%another SIZE CHECK

if size(A,1) ~= length(b)
  error('A and b must be matrices with the same number of rows.')
end

if size(P,1) ~= length(x0)
  error('P and x0 must be matrices with the same number of rows.')
end

%% initialization

iter = 0;

n = size(P,1); % number of variables

W = []; % initial set of active constraints
Wb = [];
Delta = P\(-P*x0-q); % Delta for unconstrained problem
% Delta_check = [P W';W zeros(size(W,1),size(W',2))]\[-P*x0-q; zeros(size(W,1),1)];
Lambda = [];
x = x0;

while true % iterations
    
if round(Delta,12) == 0
    
    minLambda = min(Lambda);
    if (isempty(minLambda)) || (minLambda >= 0)
        % x
        return % break while cycle
    else
        i = find(Lambda == minLambda); 
        A = [A; W(i,:)];
        b = [b; Wb(i,:)];
        W(i,:) = [];
        Wb(i,:) = [];
    end
    
else
for j = 1:size(A,1)
    checkA(j) = A(j,:)*(x+Delta) <= b(j);
end

if size(W,1) >= 1
  for j = 1:size(W,1)
    checkW(j) = W(j,:)*(x+Delta) <= Wb(j);
  end  
else
    checkW = [];
end

check = [checkA checkW];

if check == true
    Beta = 1;
    x = x + Beta*Delta;
else

    for j = 1:size(A,1)
        beta(j) = (b(j) - A(j,:)*x)/(A(j,:)*Delta);
    end

    Beta = min(beta(beta>1e-12));
    i = find(beta == Beta);
    W = [W; A(i,:)];
    A(i,:) = [];
    Wb = [Wb; b(i)];
    b(i) = [];
    x = x + Beta*Delta;

end

end

Sol = [P W';W zeros(size(W,1),size(W',2))]\[-P*x-q; zeros(size(W,1),1)];
Delta = Sol(1:n,:);
Lambda = Sol(n+1:end,:);

iter = iter+1;
checkA = []; checkW = [];

end 

end
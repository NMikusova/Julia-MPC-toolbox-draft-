%% active set test

clear
close all
clc

%% QP problem

P = [1 0;0 2];
q = [-4; -4];

N = 20; % number of edges
t = 0:2*pi/N:2*pi;
V = [(cos(t)+1)', (sin(t)+1)'];

PP = Polyhedron('V', V);
V = PP.V;

A = PP.H(:, 1:2);
b = PP.H(:, 3);
x0   = [1; 1];

[QP,J,status] = quadprog(P,q,A,b);
if status == 1
    sprintf('vsetko OK')
else
    sprintf('uz quadprog to nevie riesit, tak naco to budes robit ty')
end

%% iteration ZERO

iter = 0;

n = size(P,1); % number of variables

W = []; Wb = []; % initial set of active constraints
Delta = P\(-P*x0-q) % Delta for unconstrained problem
Delta_check = [P W';W zeros(size(W,1),size(W',2))]\[-P*x0-q; zeros(size(W,1),1)];
Lambda = [];
x = x0
% Delta = 0

%% more iterations

search = true;

while search
    
    
if round(Delta,9) == 0
    minLambda = min(Lambda)
    if (isempty(minLambda)) || (minLambda >= 0)
        x
        return
    else
        i = find(Lambda == minLambda);
        A = [A; W(i,:)];
        b = [b; Wb(i,:)];
        W(i,:) = [];
        Wb(i,:) = [];
    end
else
%% ELSE %%    
for j = 1:size(A,1)
    checkA(j) = A(j,:)*(x+Delta) <= b(j);
end
checkA

if size(W,1) >= 1
  for j = 1:size(W,1)
    checkW(j) = W(j,:)*(x+Delta) <= Wb(j);
  end  
else
    checkW = [];
end

check = [checkA checkW]

if check == true
    Beta = 1
    x = x + Beta*Delta
else

for j = 1:size(A,1)
    beta(j) = (b(j) - A(j,:)*x)/(A(j,:)*Delta);
end
beta
Beta = min(beta(beta>1e-9))
i = find(beta == Beta)
W = [W; A(i,:)]
A(i,:) = [];
A
Wb = [Wb; b(i)]
b(i) = [];
b
x = x + Beta*Delta

end

end

Sol = [P W';W zeros(size(W,1),size(W',2))]\[-P*x-q; zeros(size(W,1),1)]
Delta = Sol(1:n,:)
Lambda = Sol(n+1:end,:)

iter = iter+1
checkA = []; checkW = [];

end 



function activeset(P,q,A,b,x0)

##  Assumptions

if size(P,1) != size(P,2)
  error("P must be square matrix")
end

if minimum((eig(P))[1]) <= 0
  error("P must be positive definite matrix")
end

## SIZE CHECK

if size(P,1) != length(q)
  error("P and q must be matrices with the same number of rows.")
end

## Unconstrained Quadratic Programming

if isempty(A) && isempty(b) && isempty(x0)
  x = -inv(P)*q
  return x
end

## another SIZE CHECK

if size(A,1) != length(b)
  error("A and b must be matrices with the same number of rows.")
end

if size(P,1) != length(x0)
  error("P and x0 must be matrices with the same number of rows.")
end

## iteration ZERO

iter = 0;
n = size(P,1); # number of variables

W = Array{Float64}(0,size(A,2))
Wb = []
Delta = inv(P)*(-P*x0-q)
Lambda = []
x = x0

## iterations

while true
  if maximum(Delta) >= 0 && maximum(Delta) <= 10e-9
    if (isempty(Lambda)) || (minimum(Lambda) >= 0)
      return x, iter
    else
      i = findfirst(Lambda,minLambda)
      A = [A; W[i,:]']
      b = [b; Wb[i,:]']
      W = W[1:end .!= i, :]
      Wb = Wb[1:end .!= i, :]
    end
  else
    checkA = zeros(Int,size(A,1),1)
    checkW = zeros(Int,size(W,1),1)
    for j = 1:size(A,1)
      (A[j,:]'*(x+Delta))[1] <= b[j] ? checkA[j] = true : checkA[j] = false
    end
    if size(W,1) >= 1
      for j = 1:size(W,1)
        (W[j,:]'*(x+Delta))[1] <= Wb[j] ? checkW[j] = true : checkW[j] = false
      end
      check = vcat(checkA,checkW)
    else
      check = checkA
    end

if count(x->x==true, check) == length(check)
    Beta = 1
    x = x + Beta*Delta
else
beta_v = zeros(size(A,1),1)
for j = 1:size(A,1)
  beta_v[j] = (b[j] - A[j,:]'*x)[1]/(A[j,:]'*Delta)[1]
end
Beta = minimum(beta_v[beta_v .> 1e-9])
i = findfirst(beta_v, Beta)
W = vcat(W, A[i,:]')
A = A[1:end .!= i, :]
Wb = vcat(Wb, b[i])
b = b[1:end .!= i, :]
x = x + Beta*Delta

end
end

Sol = inv([P W';W zeros(size(W,1),size(W',2))])*[-P*x-q; zeros(size(W,1),1)]
Delta = Sol[1:n,:]
Lambda = Sol[n+1:end,:]

iter = iter+1
end
end

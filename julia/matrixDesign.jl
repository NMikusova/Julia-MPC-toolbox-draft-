function denseMPC(A,B,C,D,x0,um1,Qy,Qu,ref,N,umin,umax;ymin=[],ymax=[],dumin=[],dumax=[])

## SIZE CHECK
if size(A,1) != size(A,2)
  error("A must be square matrix")
end

if size(A,1) != size(B,1)
  error("A and B must be matrices with the same number of rows.")
end

if size(A,2) != size(C,2)
  error("A and C must be matrices with the same number of columns.")
end

if size(B,2) != size(D,2)
  error("B and D must be matrices with the same number of columns.")
end

if size(C,1) != size(D,1)
  error("C and D must be matrices with the same number of rows.")
end

## NUMBER OF STATES
nx = size(A,1)
## NUMBER OF INPUTS
nu = size(B,2)
## NUMBER OF OUTPUTS
ny = size(C,1)

## Y = Ct*x0 + Dt*U
Ct = C
for i = 1:(N-1)
  Ct = [Ct; C*(A^i)]
end

c = D
for i = 0:(N-2)
  c = [c; C*(A^i)*B]
end
Dt = c
for i = 1:(N-1)
  ci = [zeros(ny*i,nu);c[1:(end-ny*i),:]]
  Dt = [Dt ci]
end
###################################################

###################################################
## dU = L*U + lambda*u-1
lambda = [-eye(nu,nu);zeros(nu*(N-1),nu)]

c = [eye(nu,nu);-eye(nu,nu);zeros(nu*(N-2),nu)]
L = c
for i = 1:(N-1)
  if i == (N-1)
    L = [L [zeros(nu*(N-1),nu);eye(nu,nu)]]
  else
    L = [L circshift(c, i*nu)]
  end
end
###################################################

###################################################
## min (Y-R)'*Qyt*(Y-R) + dU'*Qut*dU
Qyt = kron(eye(N),Qy)
Qut = kron(eye(N),Qu)

R = []
if length(ref) == ny
  for i = 1:N
    R = [R; ref]
  end
else
  error("Reference has the wrong size")
end

###################################################

if length(umax) == nu
  Umax = []
  for i = 1:N
    Umax = [Umax; umax]
  end
  bieq = Umax
  Aieq = eye(nu*N)
else
  error("umax has the wrong size.")
end

if length(umin) == nu
  Umin = []
  for i = 1:N
    Umin = [Umin; umin]
  end
  bieq = vcat(bieq, -Umin)
  Aieq = vcat(Aieq, -eye(nu*N))
else
  error("umin has the wrong size")
end

if isempty(ymax) == false
  if length(ymax) == ny
    Ymax = []
    for i = 1:N
      Ymax = [Ymax; ymax]
    end
  else
    error("ymax has the wrong size")
  end
  bieq = vcat(bieq, Ymax-Ct*x0)
  Aieq = vcat(Aieq, Dt)
end

if isempty(ymin) == false
  if length(ymin) == ny
    Ymin = []
    for i = 1:N
      Ymin = [Ymin; ymin]
    end
  else
    error("ymin has the wrong size")
  end
  bieq = vcat(bieq, -Ymin+Ct*x0)
  Aieq = vcat(Aieq, -Dt)
end

if isempty(dumax) == false
  if length(dumax) == nu
    dUmax = []
    for i = 1:N
      dUmax = [dUmax; dumax]
    end
  else
    error("dumax has the wrong size.")
  end
  bieq = vcat(bieq, dUmax - lambda*um1)
  Aieq = vcat(Aieq, L)
end

if isempty(dumin) == false
  if length(dumin) == nu
    dUmin = []
    for i = 1:N
      dUmin = [dUmin; dumin]
    end
  else
    error("dumin has the wrong size")
  end
  bieq = vcat(bieq, -dUmin + lambda*um1)
  Aieq = vcat(Aieq, -L)
end

## P
P = Dt'*Qyt*Dt + L'*Qut*L

## qT
q = 2*(x0'*Ct'*Qyt*Dt - R'*Qyt*Dt + um1'*lambda'*Qut*L)'

## r
r = x0'*Ct'*Qyt*Ct*x0 + R'*Qyt*R + um1'*lambda'*Qut*lambda*um1 - 2*x0'*Ct'*Qyt*R

U0 = ones(N*nu,1)*((umin+umax)/2)

return P,q,r,Aieq,bieq,U0

end

##

workspace()
include("matrixDesign.jl")
include("activeset.jl")

## state-space model

A = [-0.0664 -0.4342;0.2895 0.6574]
B = [0.2895;0.2284]
C = [-2 1.5]
D = 0

x0 = [0;0];
um1 = 0;
N = 10;
Qy = 1e3;
Qu = 1e-2;
ref = 2;
umin = -2.2;
umax = 2.2;

# M = denseMPC(A,B,C,D,x0,um1,Qy,Qu,ref,N,umin,umax)
# P = 2*M[1];q = M[2];r = M[3];
# Ai = M[4];bi = M[5];
# U0 = M[6];
# iter = activeset(P,q,Ai,bi,U0)

tic()
t = 20
uu = zeros(t,1)
yy = zeros(t,1)

for k = 1:t
## QP matrices
M = denseMPC(A,B,C,D,x0,um1,Qy,Qu,ref,N,umin,umax)
P = 2*M[1];q = M[2];r = M[3];
Ai = M[4];bi = M[5];
U0 = M[6];
## activeset
iter = activeset(P,q,Ai,bi,U0)
um1 = (iter[1])[1]
uu[k] = um1
yy[k] = (C*x0)[1] + D*um1
x0 = A*x0 + B*um1
end
toc()

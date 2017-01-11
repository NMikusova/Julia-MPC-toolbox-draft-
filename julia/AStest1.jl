## Active-set test no. 1

workspace()
include("activeset.jl")

## QP problem

P = [2 0;0 2];
q = [-4;-4];

A = [2 1;1 -1;-1 -1;-2 1];
b = [2;1;1;2];
x0 = [-1;0];

## active-set

tic()
xiter = activeset(P,q,A,b,x0)
toc()

## Active-set test no. 2

workspace()
include("activeset.jl")

## QP problem

P = [1 0;0 2];
q = [-2; -6];

A = [1 1; -1 2; 2 1; -1 0; 0 -1];
b = [2; 2; 3; 0; 0];
x0 = [0.8;0.8];

## active-set

xiter = activeset(P,q,A,b,x0)

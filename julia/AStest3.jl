## Active-set test no. 2

workspace()
include("activeset.jl")

## QP problem

P = [1 0;0 2];
q = [-2; -6];

A = [];
b = [];
x0 = [];

## active-set

xiter = activeset(P,q,A,b,x0)

%% active set test

clear
close all
clc

%% QP problem

P = [2 0;0 2];
q = [-4;-4];

A = [2 1;1 -1;-1 -1;-2 1];
b = [2;1;1;2];
x0 = [-1;0];

tic;
[QP,J,status] = quadprog(P,q,A,b);
toc; QP

if status == 1
    sprintf('vsetko OK')
else
    sprintf('uz quadprog to nevie riesit, tak naco to budes robit ty')
end

tic;
[x,W,iter] = activeSet(P,q,A,b,x0);
toc; x
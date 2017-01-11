%% active set test

clear
close all
clc

%% QP problem

P = [1 0;0 2]; 
q = [-2; -6];
A = [1 1; -1 2; 2 1; -1 0; 0 -1];
b = [2; 2; 3; 0; 0];
x0 = [0.8;0.8];

[QP,J,status] = quadprog(P,q,A,b);
QP
if status == 1
    sprintf('vsetko OK')
else
    sprintf('uz quadprog to nevie riesit, tak naco to budes robit ty')
end

tic;
[x,W,iter] = activeSet(P,q,A,b,x0);
toc; x
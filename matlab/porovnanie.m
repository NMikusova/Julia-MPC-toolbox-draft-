% 
clear
close all
warning off

%%

Ts = 0.5; % perioda vzorkovania

A = [-0.0664 -0.4342;0.2895 0.6574];
B = [0.2895;0.2284];
C = [-2 1.5];
D = 0;

%% MPC

N = 10;
Qy = 1e3;
Qu = 1e-2;

%% ohranicenia
xx = cell(N,1);
uu = cell(N,1);
yy = cell(N,1);
rr = cell(N,1);

nx = 2;
nu = 1;
ny = 1;

for k = 1:N;
    xx{k} = sdpvar(nx, 1);
    uu{k} = sdpvar(nu, 1);
    yy{k} = sdpvar(ny, 1);
    rr{k} = sdpvar(ny, 1);
end
xx{k+1} = sdpvar(nx, 1);

cst =[]; 
obj = 0;                                 % consttrains
um = sdpvar(nu, 1);
for k = 1:N;
    if k == 1;
        du = uu{k} -um;
    else
        du = uu{k} - uu{k-1};
    end
    cst = cst + [xx{k+1} == A*xx{k} + B*uu{k}]; 
    cst = cst + [yy{k} == C*xx{k} + D*uu{k}];
    cst = cst + [-2.2 <= uu{k} <= 2.2];
%     cst = cst + [-0.75 <= du <= 1];
    
    obj = obj + (yy{k} - rr{k})'*Qy*(yy{k} - rr{k}) +...
        du'*Qu*du;
    
end

% options = sdpsettings('verbose',0,'solver','gurobi','cachesolver'1);

%% build reference

R = cell(40,1);
Rplot = ones(20,1) * 2;
for i = 1:length(R)
    R{i} = 2;
end

tic;
kf = 20;
x = zeros(nx, kf);
y = zeros(ny, kf);
u = zeros(nu, kf);
u_in = 0;

for k = 1:1:kf
    %MPC
    dR = cell2mat(R(k:k+N-1));
    cst_run = [];
    cst_run = cst + [xx{1} == x(:, k); um == u_in; [rr{1:end}]' == dR]; 
    sol = optimize(cst_run, obj);
    u(:, k) = value(uu{1});
    % ss model
    u_in = u(:, k);
    x(:, k+1) = A*x(:, k) + B*u(:, k);
    y(:, k)   = C*x(:, k) + D*u(:, k);
    
end
t1 = toc;
%%

time = 0:kf-1;

figure
% plot(time,y,'g')
stairs(time, Rplot,'r')
hold on
stairs(time, y,'b')
legend('skok','disketny y')
hold off, figure
umin = -2.2 + zeros(1,length(time));
umax = 2.2 + zeros(1,length(time));
plot(time,umin,time,umax)
grid on
hold on
stairs(time,u)
xlabel('t[s]')
ylabel('u')

%%

%%

A = [-0.0664 -0.4342;0.2895 0.6574];
B = [0.2895;0.2284];
C = [-2 1.5];
D = 0;

x0 = [0;0];
um1 = 0;
N = 10;
Qy = 1e3;
Qu = 1e-2;
ref = 2;
umin = -2.2; 
umax = 2.2;

% [P,q,r,Aieq,bieq,U0] = denseMPC(A,B,C,D,x0,um1,Qy,Qu,ref,N,umin,umax,[],[],[],[])
% tic;
% [QP,J,status] = quadprog(2*P,q,Aieq,bieq);
% toc; QP

tic;
xx2 = cell(101,1);
yy2 = cell(100,1);
uu2 = cell(100,1);

xx2{1} = x0;

for k = 1:20
    [P,q,r,Aieq,bieq,U0] = denseMPC(A,B,C,D,xx2{k},um1,Qy,Qu,ref,N,umin,umax,[],[],[],[]);
    QP = quadprog(2*P,q,Aieq,bieq);
    um1 = QP(1);
    uu2{k} = QP(1); % u0
    xx2{k+1} = A*xx2{k} + B*uu2{k};
    yy2{k} = C*xx2{k} + D*uu2{k};
end
t2 = toc;

time = 0:19;
up = cell2mat(uu2);
yp = cell2mat(yy2);
rp = ones(20,1)*ref;

figure
% plot(time,y,'g')
stairs(time, rp,'--')
grid on
hold on
xlabel('t[s]')
ylabel('y')
stairs(time, yp,'b')
legend('reference','output')
hold off, 
figure
umax = 2.2 + zeros(1,length(time));
umin = -2.2 + zeros(1,length(time));
plot(time,umax)
hold on
stairs(time,up)
grid on
plot(time,umin)
xlabel('t[s]')
ylabel('u')
legend('constraints','control action')
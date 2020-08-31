close all
clear all
clc

import casadi.*
addpath('../OPTCON')

%% Setup the general parameters
dt=0.01;
horizon=1;
T=horizon/dt;

% Setup the parameter for the robot arm
m1=1; %mass of the 1st link
m2=1; %mass of the 2nd link
l1=1; %length of the 1st link
l2=1; %length of the 2nd link
lc1=0.5; %CoM of the 1st link
lc2=0.5; %CoM of the 2nd link
I1=0.5; %moment of initia of the 1st link
I2=0.5; %moment of initia of the 2nd link
g=9.8; %gravity costant


% Declare model variables for the two-link robot arm
q1 = SX.sym('q1');
q2 = SX.sym('q2');
dq1= SX.sym('dq1');
dq2= SX.sym('dq2');
x = [q1;q2;dq1;dq2];
u1=SX.sym('tau1');
u2=SX.sym('tau2');
u = [u1;u2];
n=length(x);
m=length(u);

%% setup the system dynamics model
a1=m1*lc1^2+m2*(l1*l1+lc2*lc2)+I1+I2;
a2=m2*l1*lc2;
a3=m2*(lc2*lc2)+I2;
b1=(m1*lc1+m2*l1);
b2=m2*lc2;

M11=a1+2*a2*cos(q2);
M12=a3+a2*cos(q2);
M21=M12;
M22=a3;
invM11=M22/(M11*M22-M12*M21);
invM12=-M12/(M11*M22-M12*M21);
invM21=-M21/(M11*M22-M12*M21);
invM22=M11/(M11*M22-M12*M21);
C11=-a2*dq2*sin(q2);
C12=-a2*(dq1+dq2)*sin(q2);
C21=a2*dq1*sin(q2);
C22=0;
C=[C11,C12;C21,C22];
G1=b1*g*cos(q1)+b2*g*cos(q1+q2);
G2=b2*g*cos(q1+q2);
G=[G1;G2];
ddq1=[invM11, invM12]*(-C*[dq1;dq2]-G+u);
ddq2=[invM21, invM22]*(-C*[dq1;dq2]-G+u);
f = [q1+dq1*dt;
     q2+dq2*dt;
     dq1+dt*ddq1;
     dq2+dt*ddq2];
% Discrete time dynamics
dyn = Function('dyn', {x, u}, {f}, {'X','U'}, {'f'});
dfxFun = Function('dfx',{x,u},{jacobian(f,x)});
dfuFun = Function('dfu',{x,u},{jacobian(f,u)});

%% setup the cost function
features = [u1^2;
            u2^2;
            u1*u2];
weights=[0.6,0.3, 0.1]';
% cost function
cost= Function('cost',{x,u},{weights'*features}, {'X','U'}, {'c'});


%% use the oc solver the solve the optimal control probme
x0=[0,0,0,0]';
xgoal=[pi/2,-pi/2,0,0]';
sol=OCsolver_FixEnd(x0,xgoal,T,dyn,cost);
clc


%% do the inverse optimal control with 3 features
features_1=[u1^2, u2^2, u1*u2]';
r=3;
% cost function
dpxFun_1=Function('dpx',{x,u},{jacobian(features_1,x)});
dpuFun_1=Function('dpu',{x,u},{jacobian(features_1,u)});

st=50;
l=1;
[dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun_1,dpuFun_1,sol);
H2=dfu'*dfx';
H1=dfu'*dpx'+dpu';
H=[H1,H2];
rankH_1=rank(H);
% solve the weights
results_1=SolveH(H,r);
error_1=SolveError(SolveH(H,r),[0.6 0.3 0.1]);
for l=2:50
    [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun_1,dpuFun_1,sol);
    H1=[H1+H2*dpx';
        dfu'*dpx'+dpu'];
    H2=[H2*dfx';
        dfu'*dfx'];
    H=[H1 H2];
    % solve the weights
    error_1(end+1)=SolveError(SolveH(H,r),[0.6 0.3 0.1]);
    rankH_1(end+1)=rank(H);
end

%% do the inverse optimal control with 4 features
features_2=[u1^2, u2^2, u1*u2, u1^3, u2^3, u1*u1*u2]';
r=6;
% cost function
dpxFun_2=Function('dpx',{x,u},{jacobian(features_2,x)});
dpuFun_2=Function('dpu',{x,u},{jacobian(features_2,u)});

st=50;
l=1;
[dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun_2,dpuFun_2,sol);
H2=dfu'*dfx';
H1=dfu'*dpx'+dpu';
H=[H1,H2];
rankH_2=rank(H);
% solve the weights
error_2=SolveError(SolveH(H,r),[0.6 0.3 0.1, 0,0, 0]);
for l=2:50
    [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun_2,dpuFun_2,sol);
    H1=[H1+H2*dpx';
        dfu'*dpx'+dpu'];
    H2=[H2*dfx';
        dfu'*dfx'];
    H=[H1 H2];
    % solve the weights
    error_2(end+1)=SolveError(SolveH(H,r),[0.6 0.3 0.1, 0,0, 0]);
    rankH_2(end+1)=rank(H);
end


%% do the inverse optimal control with 6 features
features_3=[u1^2, u2^2, u1*u2, u1^3, u2^3,  u1*u2*u2, u1^4, u2^4, u2^2*u1^2, u1^3*u2]';
r=10;
% cost function
dpxFun_3=Function('dpx',{x,u},{jacobian(features_3,x)});
dpuFun_3=Function('dpu',{x,u},{jacobian(features_3,u)});

st=50;
l=1;
[dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun_3,dpuFun_3,sol);
H2=dfu'*dfx';
H1=dfu'*dpx'+dpu';
H=[H1,H2];
rankH_3=rank(H);
% solve the weights and error
error_3=SolveError(SolveH(H,r),[0.6 0.3 0.1, 0,0, 0, 0, 0, 0, 0]);
for l=2:50
    [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun_3,dpuFun_3,sol);
    H1=[H1+H2*dpx';
        dfu'*dpx'+dpu'];
    H2=[H2*dfx';
        dfu'*dfx'];
    H=[H1 H2];
    % solve the weights
    error_3(:,end+1)=SolveError(SolveH(H,r),[0.6 0.3 0.1, 0,0, 0, 0, 0, 0, 0]);
    rankH_3(end+1)=rank(H);
end


%% do the plot
figure(1)
subplot(2,1,1)
plot(1:length(rankH_1),rankH_1,'LineWidth',3,'LineStyle','-')
hold on
plot(1:length(rankH_2),rankH_2,'LineWidth',3,'LineStyle','--')
plot(1:length(rankH_3),rankH_3,'LineWidth',3,'LineStyle','-.')
legend('$\{\tau_1^2,\tau_2^2, \tau_1\tau_2\}$',...
    '$\{\tau_1^2,\tau_2^2, \tau_1\tau_2, \tau_1^3,\tau_2^3,  \tau_1^2\tau_2\}$',...
    '$\{\tau_1^2,\tau_2^2, \tau_1\tau_2, \tau_1^3,\tau_2^3,  \tau_1\tau_2^2, \tau_1^4,\tau_2^4, \tau_1^3\tau_2, \tau_1^2\tau_2^2  \}$','Interpreter','latex')
xlim([0 50])
ylim([1 14])
ylabel('rank H')
box on
grid on
subplot(2,1,2)
plot(1:length(rankH_1),error_1,'LineWidth',3,'LineStyle','-')
hold on
plot(1:length(rankH_2),error_2,'LineWidth',3,'LineStyle','--')
plot(1:length(rankH_3),error_3,'LineWidth',3,'LineStyle','-.')
legend('$\{\tau_1^2,\tau_2^2, \tau_1\tau_2\}$',...
    '$\{\tau_1^2,\tau_2^2, \tau_1\tau_2, \tau_1^3,\tau_2^3,  \tau_1^2\tau_2\}$',...
    '$\{\tau_1^2,\tau_2^2, \tau_1\tau_2, \tau_1^3,\tau_2^3,  \tau_1\tau_2^2, \tau_1^4,\tau_2^4, \tau_1^3\tau_2, \tau_1^2\tau_2^2 \}$','Interpreter','latex')
ylabel('$e_{\omega}$', 'Interpreter','latex')
xlabel('observation length $l$ (with starting time $t=50$)','Interpreter','latex')
xlim([0 50])
box on
grid on




%% solve the differentiable dynamics and cost function
function [dfx,dfu,dpx,dpu]=DiffDynCost(t,dfxFun,dfuFun,dpxFun,dpuFun,sol)
% note here index for u is t-1, index for x is t
dfu=full(dfuFun(sol.x(:,t),sol.u(:,t)));
dpu=full(dpuFun(sol.x(:,t),sol.u(:,t)));
dfx=full(dfxFun(sol.x(:,t+1),sol.u(:,t+1)));
dpx=full(dpxFun(sol.x(:,t+1),sol.u(:,t+1)));
end

%% solve the recovery matrix
function x=SolveH(H,r)
[U,D,V]=svd(H);
v=V(:,end);
v=v(1:r);
v=v/sign(v(2));
x=v/sum(v);
% options = optimoptions('quadprog','Display','off');
% n=size(H,2);
% Aeq=ones(1,n);
% Aeq(r+1:end)=0;
% beq=1;
% lb=zeros(1,n);
% lb(r+1:end)=-inf;
% ub=ones(1,n);
% ub(r+1:end)=inf;
% x=quadprog(H'*H,[],[],[],Aeq,beq,lb,ub,[],options);
% x=x(1:r);
end

%% evaluate the error index
function e=SolveError(v,vg)
if size(v,2)>1; v=v';end
if size(vg,2)>1; vg=vg'; end
c=dot(v,vg)/dot(v,v);
e=norm(c*v-vg)/norm(vg);

end



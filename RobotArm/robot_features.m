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

% add the noise level
sigma=1e-4;
sol.x=sol.x+sigma*randn(size(sol.x));
sol.u=sol.u+sigma*randn(size(sol.u));


%% do the inverse optimal control with 3 features
features=[u1^2, u2^2, u1*u2,u1^3, u2^3, u1*u1*u2, u1*u2*u2, ]';
r=length(features);
true_weights=zeros(r,1);
true_weights(1:3)=weights;
% cost function
dpxFun=Function('dpx',{x,u},{jacobian(features,x)});
dpuFun=Function('dpu',{x,u},{jacobian(features,u)});

gamma=45;

st=-1;
results=[];
stime=[];
oblength=[];
while st<=T
    st=st+1;
    l=1;
    [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun,dpuFun,sol);
    H2=dfu'*dfx';
    H1=dfu'*dpx'+dpu';
    H=[H1,H2];
    kappa=EvaluateRank(H);
    l=l+1;  % start from 2
    while kappa<gamma && st+l<=T
        [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun,dpuFun,sol);
        H1=[H1+H2*dpx';
            dfu'*dpx'+dpu'];
        H2=[H2*dfx';
            dfu'*dfx'];
        H=[H1 H2];
        kappa=EvaluateRank(H);
        l=l+1;
    end
    
    if (st+l)<=T
        stime(end+1)=st;
        results(:,end+1)=SolveH(H,r);
        oblength(end+1)=l;
        clc % for display
        current_time=st
        current_length=l
        current_result=SolveH(H,r)
    else
        break
    end
    
    
    
end

%% results data post-processing
error=[];
for k=1:length(stime)
    error(end+1)=SolveError(results(:,k),true_weights);
end
gamma
meanerror=mean(error)
meanlength=mean(oblength)/100





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

%% evaluate the rank index 
function gamma=EvaluateRank(H)
[r,c]=size(H);
if r<2*c
    gamma=0;
else
    s=svd(H);
    s=sort(s);
    gamma=s(2)/s(1);
end
end




close all
clear all
clc
import casadi.*
addpath('../OPTCON')

%% system setting
A=[-1 1;
  	0 1];
B=[1,3]';
% define the state and input variable
x=SX.sym('x',2);
u=SX.sym('u',1);
% define the dynamics
f=A*x+B*u;
dyn = Function('dyn', {x, u}, {f}, {'X','U'}, {'f'});
dfxFun = Function('dfx',{x,u},{jacobian(f,x)});
dfuFun = Function('dfu',{x,u},{jacobian(f,u)});

%% setup the cost function
features=[x(1)^2, x(2)^2, u^2]';
weights=[0.1,0.3,0.6]';
% cost function
cost= Function('cost',{x,u},{weights'*features}, {'X','U'}, {'c'});

%% use the oc solver the solve the optimal control probme
x0=[2,-2]';
T=50;
sol=OCsolver_FreeEnd(x0,T,dyn,cost);
clc

%% do the inverse optimal control
features=[x(1)^2, x(2)^2, u^2]';
r=3;
% cost function
phi=Function('feature',{x, u}, {features}, {'X','U'}, {'phi'});
dpxFun=Function('dpx',{x,u},{jacobian(features,x)});
dpuFun=Function('dpu',{x,u},{jacobian(features,u)});

st=2;
l=1;
[dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun,dpuFun,sol);
H2=dfu'*dfx';
H1=dfu'*dpx'+dpu';
H=[H1,H2];
rankH=rank(H);
% solve the weights
results=SolveH(H,r);
restuls_inverseKKT=SolveH(H1,r);
for l=2:T-st
    [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun,dpuFun,sol);
    H1=[H1+H2*dpx';
        dfu'*dpx'+dpu'];
    H2=[H2*dfx';
        dfu'*dfx'];
    H=[H1 H2];
    H=0.001*H;
%     H=H/norm(H,'fro');
    % solve the weights
    results(:,end+1)=SolveH(H,r);
    restuls_inverseKKT(:,end+1)=SolveH(H1,r);
    rankH(end+1)=rank(H);
end

error=Compute_Error(results',weights');
error_inverseKKT=Compute_Error(restuls_inverseKKT',weights');



figure(1)
plot(1:length(error),error,'LineWidth',3)
hold on
plot(1:length(error_inverseKKT),error_inverseKKT,'LineWidth',3)
legend('Proposed method', 'Inverse-KKT method')
ylabel('$e_{\omega}$','Interpreter','latex')
grid on
box on
ylim([0,1])
xlabel('Observation length $l$ ($t=40$)','FontWeight','bold','Interpreter','latex')













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


%% compute the error
function error=Compute_Error(estcosts,truecost)

[l,n]=size(estcosts);
cost=truecost/norm(truecost);

error=[];
for i=1:l
    est=estcosts(i,:);
    c=dot(est,cost)/dot(est,est);
    e=norm(c*est-cost);
    error(end+1)=e;
end

end



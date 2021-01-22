close all
clear all
clc
import casadi.*
addpath('../OPTCON')

%% system setting
A=[-1 1;
  	0 1];
B=[1,3]';
Q=[0.1, 0;
    0, 0.3];
R=0.6;
x0=[2,-2]';

% solve the LQR problem
[K,S,e]=dlqr(A,B,Q,R);
K

% solve the trajectory
traj_x=x0;
traj_u=[];


for k=1:1000
    traj_u(:,end+1)=-K*traj_x(:,end);
    traj_x(:,end+1)=A*traj_x(:,end)+B*traj_u(:,end);
end

sol.x=traj_x;
sol.u=traj_u;

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



%% do the inverse optimal control
features=[x(1)^2, x(2)^2, u^2]';
r=3;
% cost function
phi=Function('feature',{x, u}, {features}, {'X','U'}, {'phi'});
dpxFun=Function('dpx',{x,u},{jacobian(features,x)});
dpuFun=Function('dpu',{x,u},{jacobian(features,u)});

st=8;
l=1;
[dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun,dpuFun,sol);
H2=dfu'*dfx';
H1=dfu'*dpx'+dpu';
H=[H1,H2];
rankH=rank(H);
% solve the weights
results=SolveH(H,r);
for l=2:25
    [dfx,dfu,dpx,dpu]=DiffDynCost(st+l,dfxFun,dfuFun,dpxFun,dpuFun,sol);
    H1=[H1+H2*dpx';
        dfu'*dpx'+dpu'];
    H2=[H2*dfx';
        dfu'*dfx'];
    H=[H1 H2];
    % solve the weights
    results(:,end+1)=SolveH(H,r);
    rankH(end+1)=rank(H);
end

figure(1)
subplot(2,1,1)
plot(1:length(rankH),rankH,'LineWidth',3)
ylabel('Rank H','FontWeight','bold','Interpreter','latex')
grid on
box on
xlim([0 length(rankH)])
ylim([1,5])

results(1,results(1,:)>1)=0.8;
results(2,results(2,:)<-0.2)=0;

subplot(2,1,2)
plot(1:length(rankH),results(1,:),'LineWidth',3,'LineStyle','-')
hold on
plot(1:length(rankH),results(2,:),'LineWidth',3,'LineStyle',':')
plot(1:length(rankH),results(3,:),'LineWidth',3,'LineStyle','--')
legend('$\hat\omega_1$','$\hat\omega_2$','$\hat\omega_3$','Interpreter','latex')
ylabel('$\hat\omega$','FontWeight','bold','Interpreter','latex')
grid on
box on
xlim([0 length(rankH)])
ylim([-0.2,0.8])
yticks([-0.2 0 0.2 0.4 0.6 0.8])
xlabel('Observation length $l$ (with starting time $t=8$)','FontWeight','bold','Interpreter','latex')





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



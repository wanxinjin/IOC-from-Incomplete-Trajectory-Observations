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
weights=[0.6,0.3,0.1]';
% cost function
phi=Function('feature',{x, u}, {features}, {'X','U'}, {'phi'});
cost= Function('cost',{x,u},{weights'*features}, {'X','U'}, {'c'});
dpxFun=Function('dpx',{x,u},{jacobian(features,x)});
dpuFun=Function('dpu',{x,u},{jacobian(features,u)});

%% use the oc solver the solve the optimal control probme
x0=[2,-2]';
T=50;
sol=OCsolver_FreeEnd(x0,T,dyn,cost);
traj_x=sol.x;
traj_u=sol.u;
traj_t=0:T+1;


%% do the plot 
figure(1)
subplot(2,1,1)
plot(0:T+1,traj_x(1,:),'LineWidth',3)
xlim([0,51])
hold on
plot(0:T+1,traj_x(2,:),'LineWidth',3)
grid on
ylabel('$x$','interpreter','latex')
legend('$x_1$', '$x_2$','interpreter','latex')
subplot(2,1,2)
plot(0:T,traj_u,'LineWidth',3)
xlim([0,51])
grid on
ylabel('$u$','interpreter','latex')
xlabel('time')
legend('$u$','interpreter','latex')


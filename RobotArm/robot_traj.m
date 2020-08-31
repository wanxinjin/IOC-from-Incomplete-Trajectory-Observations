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
phi = Function('feature', {x, u}, {features}, {'X','U'}, {'phi'});
cost= Function('cost',{x,u},{weights'*features}, {'X','U'}, {'c'});
dpxFun=Function('dpx',{x,u},{jacobian(features,x)});
dpuFun=Function('dpu',{x,u},{jacobian(features,u)});

%% use the oc solver the solve the optimal control probme
x0=[0,0,0,0]';
xgoal=[pi/2,-pi/2,0,0]';
sol=OCsolver_FixEnd(x0,xgoal,T,dyn,cost);
clc


%% do the plot 
figure(1)
subplot(2,1,1)
plot(0:T+1,sol.x(1,:),'LineWidth',3)
xlim([0,T+1])
hold on
plot(0:T+1,sol.x(2,:),'LineWidth',3)
plot(0:T+1,sol.x(3,:),'LineWidth',3)
plot(0:T+1,sol.x(4,:),'LineWidth',3)
grid on
ylabel('$x$','interpreter','latex')
legend('$\theta_1$', '$\theta_2$', '$\dot{\theta}_1$', '$\dot{\theta}_2$','interpreter','latex')

subplot(2,1,2)
plot(0:T,sol.u(1,:),'LineWidth',3)
hold on
plot(0:T,sol.u(2,:),'LineWidth',3)
grid on
ylabel('$u$','interpreter','latex')
xlabel('Time')
legend('$\tau_1$','$\tau_2$', 'interpreter','latex')

%190010047
%assignment-2
close all
clc 
load System1_Parameters.mat
%system parameters 
U_eq = sys.Us;
X_eq = sys.Xs;
Y_eq = sys.Ys;
Ds = sys.Ds;
alpha = sys.alfa;
dk_sigma = sys.dk_sigma;
meas_sigma = sys.meas_sigma;

%matrices for storing states and inputs 
Ns = 100;  %no of samples 

n_states = 3;
n_inputs = 2; 
n_outputs = 2;

X = zeros(n_states,Ns);  %states 
U = zeros(n_inputs,Ns);  %Inputs 
Y = zeros(n_outputs,Ns);  %Outputs 
eita = zeros(n_outputs,Ns);

x = zeros(n_states,Ns);  %perturbed states 
y = zeros(n_inputs,Ns);  %perturbed outputs 
u = zeros(n_inputs,Ns);  %perturbed inputs 
Disturbance = zeros(1,Ns); 
Refer_deviation = zeros(n_outputs,Ns);
error = zeros(n_outputs,Ns);
C_pi = [8.5396/0.7278 0;0 -1.5703/0.5286];
D_pi = [8.5396 0;0 -1.5703];
C = [0,1,0;0,0,1];
T = zeros(1,Ns);  %array for time 

%initial condition
X(:,1) = X_eq;
Y(:,1) = Y_eq;
sampl_T = 0.1;

%input Constraints 
U_low = [0 0]';
U_high = [250 60]';

%Noice Switch 
Noise_ON = 0;

Z = zeros(3,Ns);
%%%%%%%% servo control problem %%%%%%%% 
for i = 1:Ns-1
    Z(:,i) = X_eq;
    T(i) = sampl_T*(i-1);
    if i < 5
        Refer_deviation(:,i) = [0,0]';
    else 
        Refer_deviation(:,i) = [-5,0]';  
    end
    Disturbance(i) = Ds+Noise_ON*normrnd(0,0.015);
    %%%%% PI controller %%%%%% 
    y(:,i) = Y(:,i) - Y_eq;
    error(:,i) = Refer_deviation(:,i) - y(:,i);
    eita(:,i+1) = eita(:,i) + sampl_T*error(:,i);
    u(:,i) = C_pi*eita(:,i) + D_pi*error(:,i);
    U(:,i) = u(:,i) + U_eq;
    %%%%% Input Saturation %%%%%%%
    for j = 1:n_inputs
        if U(j,i)>=U_low(j) && U(j,i) <= U_high(j)
            continue
        elseif  U(j,i) < U_low(j)
            U(j,i) = U_low(j);
        elseif  U(j,i) > U_high(j)
            U(j,i) = U_high(j);
        end
    end
    %%%%% Dynamic-Simulation %%%%%%%
    n = 10;
    X0 = X(:,i);
    h = sampl_T/n;
    for m = 1:n
        X1 = X0 + h*System1_Dynamics(X0,U(:,i),Disturbance(i),alpha);
        X0 = X1;
    end 
    X(:,i+1) = X1;
    V = [normrnd(0,0.2) normrnd(0,0.25)]';
    Y(:,i+1) = C*X(:,i+1) + Noise_ON*V;

end
T(Ns) = 10;
Z(:,Ns) = X_eq;
Refer_deviation(:,Ns) = [-5,0]';
U(:,Ns) = U(:,Ns-1);
Disturbance(Ns) = Ds+Noise_ON*normrnd(0,0.015);
%%%%% Plots %%%%%%
figure(1)
subplot(3,1,1) 
plot(T,X(1,:)),grid on ,xlabel('Time in minutes '),ylabel('X1')
hold on 
plot(T,Z(1,:))
hold off
subplot(3,1,2) 
plot(T,X(2,:)),grid on ,xlabel('Time in minutes '),ylabel('X2')
hold on 
plot(T,Z(2,:))
hold off
subplot(3,1,3) 
plot(T,X(3,:)),grid on ,xlabel('Time in minutes '),ylabel('X3')
hold on 
plot(T,Z(3,:))
hold off
figure(2)
subplot(2,1,1)
plot(T,Y(1,:)),grid on ,xlabel('Time in minutes '),ylabel('Y1')
hold on 
plot(T,Refer_deviation(1,:)+Y_eq(1))
hold off
subplot(2,1,2)
plot(T,Y(2,:)),grid on ,xlabel('Time in minutes '),ylabel('Y2')
hold on 
plot(T,Refer_deviation(2,:)+Y_eq(2))
hold off
figure(3)
subplot(3,1,1)
plot(T,U(1,:)),grid on ,xlabel('Time in minutes '),ylabel('U1')
subplot(3,1,2)
plot(T,U(2,:)),grid on ,xlabel('Time in minutes '),ylabel('U2')
subplot(3,1,3)
plot(T,Disturbance),grid on ,xlabel('Time in minutes '),ylabel('Disturbance')
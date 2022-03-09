%190010047
%assignment-3
% State Augmentation Algorithm 
close all
clc
load System1_Continuous_LinMod.mat
load System1_Parameters.mat
SetGraphics
%%linear pertubed model matrices
A = A_mat;
B = B_mat;
C = C_mat;
D = zeros(2,2);
H = H_mat;
sampl_T = 0.1;

%system parameters 
U_eq = sys.Us;
X_eq = sys.Xs;
Y_eq = sys.Ys;
Ds = sys.Ds;
alpha = sys.alfa;
dk_sigma = sys.dk_sigma;
meas_sigma = sys.meas_sigma;
%%Discrete Time Linear model Pertubation matrices

syms s 
%Laplace inverse of inverse of (SI-A) calculated at sampling time 

phi = double(subs(ilaplace(inv(s*eye(3,3)-A)),sampl_T));

% when A is invertible Tau matrix is (Phi-I)*inv(A)*B

tau_u = (phi - eye(3))*(A\B);
tau_h = (phi - eye(3))*(A\H);

%matrices for storing states and inputs 
Ns = 250;
n_states = 3;
n_inputs = 2; 
n_outputs = 2;

X = zeros(n_states,Ns);  %states 
x = zeros(n_states,Ns);  %perturbed states
x_estimated_augmented  = zeros(n_states+n_inputs,Ns);  %estimated states
x_s = zeros(n_states,Ns); %steady states

U = zeros(n_inputs,Ns);  %Inputs 
u = zeros(n_inputs,Ns);  %perturbed inputs 
u_s = zeros(n_inputs,Ns);  %steady inputs 

Y = zeros(n_outputs,Ns);  %Outputs 
y = zeros(n_inputs,Ns);  %perturbed outputs 
 
Disturbance = zeros(1,Ns); 
Reference_deviation = zeros(n_outputs,Ns);
error = zeros(n_outputs,Ns);
T = zeros(1,Ns);  %array for time 

%initial condition
X(:,1) = X_eq +[0.005 2 -2]';
Y(:,1) = C*X(:,1);

%input Constraints 
U_low = [0 0]';
U_high = [250 60]';

%augmented system 
phi_augmented = [phi,tau_u;zeros(n_outputs,n_states) ,eye(n_outputs,n_inputs)];
tau_u_augmented = [tau_u ;zeros(n_outputs,n_inputs)];
C_augmented = [C ,zeros(n_outputs,n_inputs)];

%poles 
Controller_poles = [0.7 0.4 0.2]';
Observer_poles = [0.3 0.4 0.45 0.5 0.6]';

%Gain Matrices 
rank_C = rank(ctrb(phi,tau_u));
fprintf('\n Rank of Controllability Matrix: %d \n ',rank_C);
G = place(phi,tau_u,Controller_poles);
rank_O = rank(obsv(phi_augmented,C_augmented));
fprintf('\n Rank Augmented Observability Matrix: %d \n ',rank_O);
L_augmented = place(phi_augmented',C_augmented',Observer_poles);
L_augmented = L_augmented';
% Steady State Matrices 
m = (eye(n_states)-phi)\tau_u;
K_u = C*m;
for i = 1:Ns-1 
    T(i) = sampl_T*(i-1);
    if i*sampl_T >= 6 && i*sampl_T < 12
        Reference_deviation(:,i) = [-5 0]';
    end 
    if i*sampl_T >= 18 
        Disturbance(i) = 0.1;
    end 
    %%%% Control Law %%%%%%%%
    u_s(:,i) = K_u\(Reference_deviation(:,i))-x_estimated_augmented(4:end,i);
    muk = (eye(n_states)-phi)\tau_u;
    esh = (K_u)\Reference_deviation(:,i);
    x_s(:,i) = muk*esh;
    u(:,i) = u_s(:,i) - G*(x_estimated_augmented(1:3,i)-x_s(:,i));
    U(:,i) = u(:,i)+U_eq;
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
    Disturbance(i) = Disturbance(i) +Ds;
    %%%%% Dynamic-Simulation for getting measurements %%%%%%%
    n = 10;
    X0 = X(:,i);
    h = sampl_T/n;
    for m = 1:n
        X1 = X0 + h*System1_Dynamics(X0,U(:,i),Disturbance(i),alpha);
        X0 = X1;
    end 
    X(:,i+1) = X1;
    Y(:,i+1) = C*X(:,i+1);
    
    %%%%%%%%% State Augmentation Approach %%%%%%%%%%
    
    y(:,i) = Y(:,i) - Y_eq;
    error(:,i) = y(:,i) - C_augmented*x_estimated_augmented(:,i);
    
    x_estimated_augmented(:,i+1) = phi_augmented*x_estimated_augmented(:,i)+tau_u_augmented*u(:,i)+L_augmented*error(:,i); 

end 
T(Ns) =25;
U(:,Ns) = U(:,Ns-1);
x_s(:,Ns) = x_s(:,Ns-1);
u_s(:,Ns) = u_s(:,Ns-1);
error(:,Ns) = Y(:,Ns) - Y_eq- C_augmented*x_estimated_augmented(:,Ns);
Disturbance(Ns) = Ds+0.1;
%%%%%%Plots%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)
plot(T,x_s(1,:)),grid on ,xlabel('Time  '),ylabel('xs_1'),title('Target Model Steady State');
subplot(3,1,2)
plot(T,x_s(2,:)),grid on ,xlabel('Time '),ylabel('xs_2');
subplot(3,1,3)
plot(T,x_s(3,:)),grid on ,xlabel('Time  '),ylabel('xs_3');
figure(2)
subplot(2,1,1)
plot(T,u_s(1,:)),grid on ,xlabel('Time  '),ylabel('us_1'),title('Target Model Manipulated Input');
subplot(2,1,2)
plot(T,u_s(2,:)),grid on ,xlabel('Time  '),ylabel('us_2');
figure(3)
subplot(2,1,1)
plot(T,error(1,:)),grid on ,xlabel('Time '),ylabel('e_1 '),title('Measurement Error ');
subplot(2,1,2)
plot(T,error(2,:)),grid on ,xlabel('Time '),ylabel('e_2');
figure(4)
subplot(2,1,1)
stairs(T,U(1,:),"LineWidth",2),grid on ,xlabel('Time '),ylabel('U1'),title('Manipulated Inputs');
subplot(2,1,2)
stairs(T,U(2,:),"LineWidth",2),grid on ,xlabel('Time '),ylabel('U2');
figure(5)
subplot(2,1,1)
plot(T,Y(1,:)),grid on ,xlabel('Time '),ylabel('Y_1'),title('Measured Output');
hold on 
plot(T,Reference_deviation(1,:)+Y_eq(1));
hold off
subplot(2,1,2)
plot(T,Y(2,:)),grid on ,xlabel('Time '),ylabel('Y_2');
hold on 
plot(T,Reference_deviation(2,:)+Y_eq(2));
hold off
figure(6)
subplot(3,1,1)
plot(T,X(1,:)),grid on ,xlabel('Time '),ylabel('X_1'),title('State Variables');
subplot(3,1,2)
plot(T,X(2,:)),grid on ,xlabel('Time '),ylabel('X_2');
subplot(3,1,3)
plot(T,X(3,:)),grid on ,xlabel('Time '),ylabel('X_3');
figure(7)
stairs(T,Disturbance,"LineWidth",2),grid on ,xlabel('Time '),ylabel('Dk'),title('Disturbance');
figure(8)
subplot(3,1,1)
plot(T,X(1,:)-x_estimated_augmented(1,:)),grid on ,xlabel('Time  '),ylabel('eps_1'),title('State Estimation Error');
subplot(3,1,2)
plot(T,X(2,:)-x_estimated_augmented(1,:)),grid on ,xlabel('Time '),ylabel('eps_2');
subplot(3,1,3)
plot(T,X(3,:)-x_estimated_augmented(1,:)),grid on ,xlabel('Time  '),ylabel('eps_3');
figure(9)
subplot(2,1,1)
plot(T,x_estimated_augmented(4,:)),grid on ,xlabel('Time '),ylabel('beta_1'),title('Input Bias');
subplot(2,1,2)
plot(T,x_estimated_augmented(5,:)),grid on ,xlabel('Time '),ylabel('beta_2');
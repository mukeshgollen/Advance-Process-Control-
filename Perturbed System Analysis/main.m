%190010047
%assignment-1
close all
clc
load System1_Parameters.mat
load System1_Continuous_LinMod.mat
randn('seed',0);
SetGraphics
%linear pertubed model matrices
A = A_mat;
B = B_mat;
C = C_mat;
D = zeros(2,2);
H = H_mat;
%system parameters 
U_eq = sys.Us;
X_eq = sys.Xs;
Y_eq = sys.Ys;
Ds = sys.Ds;
alpha = sys.alfa;
dk_sigma = sys.dk_sigma;
meas_sigma = sys.meas_sigma;

%matrices for storing states and inputs 
Ns = 500; %no of samples 
X = zeros(3,Ns);
U = zeros(2,Ns);
X_li =zeros(3,Ns);
x = zeros(3,Ns);
Disturbance = zeros(1,Ns);
T = zeros(1,Ns);%array for time 
%initial condition
X(:,1) = X_eq;
X_li(:,1) = X_eq;
sampl_T = 0.1;
%Discrete Time Linear model Pertubation matrices

syms s 
%Laplace inverse of inverse of (SI-A) calculated at sampling time 

phi = double(subs(ilaplace(inv(s*eye(3,3)-A)),sampl_T));

% when A is invertible Tau matrix is (Phi-I)*inv(A)*B

tau_u = (phi - eye(3))*(A\B);
tau_h = (phi - eye(3))*(A\H);

%random binary input generation 
RBS_period = [25 30]';
RBS_amp = [10 4]';
m = zeros(2,1);%variable to store input 

%Dynamic simulation 
for k = 1:Ns-1
    T(k) = (k-1)*sampl_T;
    %%%% Manipulated input %%%%%%
    for i =1:2
        if (rem(k,RBS_period(i))==0)
            m(i) = RBS_amp(i)*sign(randn);
        end
    end 
    %%%% Disturbance %%%%%%%%
    if k < 250
       d = 0;
    else 
       d = -0.2;
    end
    %%%%%%%% Linear %%%%%%%%%%
    U(:,k) = m;
    Disturbance(k) = d;
    x(:,k+1) = phi*x(:,k) + tau_u*U(:,k)+tau_h*Disturbance(k);
    X_li(:,k+1) = X_eq + x(:,k+1);
    U(:,k) =  U(:,k) + U_eq;
    Disturbance(k) = Disturbance(k) + Ds;
    %%%%%%%%%%%%% Non - Linear %%%%%%%%%%%%%%%%%
 
    
    n = 10;
    X0 = X(:,k);
    h = sampl_T/n;
    for i = 1:n
        X1 = X0 + h*System1_Dynamics(X0,U(:,k),Disturbance(k),alpha);
        X0 = X1;
    end 
    X(:,k+1) = X1;
    % The above method has been used since the value of the derivative was
    % changing abruptly in the sampling interval and eventually the
    % previous method has to be updated
    
 end 
T(Ns) = 50;%time at the end of simulation 
U(:,Ns) = m + U_eq;
Disturbance(Ns) = d + Ds;

%%%%%%% Plots %%%%%%%%%%%%%%
figure(1)
plot(T,X_li(1,:)),grid on ,xlabel('Time in minutes '),ylabel('X1')
hold on 
plot(T,X(1,:)),legend('Linear','Non-Linear')
hold off 
figure(2)
plot(T,X_li(2,:)),grid on ,xlabel('Time in minutes '),ylabel('X2')
hold on 
plot(T,X(2,:)),legend('Linear','Non-Linear')
hold off 
figure(3)
plot(T,X_li(3,:)),grid  on ,xlabel('Time in minutes '),ylabel('X3')
hold on 
plot(T,X(3,:)),legend('Linear','Non-Linear')
hold off 
figure(4)
plot(T,Disturbance),grid on ,xlabel('Time in minutes '),title('Distubance')
figure(5)
subplot(2,1,1),plot(T,U(1,:)),grid on ,ylabel('U1')
subplot(2,1,2),plot(T,U(2,:)),grid on ,xlabel('Time in minutes '),ylabel('U2')




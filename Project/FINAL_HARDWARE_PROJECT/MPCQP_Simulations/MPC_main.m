
clear all 
clc 

load TCL_LQOC_AMX_24July2020.mat
load TCL_SS_Model_AMX_24July2020.mat
load TCL_Model_Parameters.mat
global TCL 

N_pred = 75;
N_con = 6;
n_states = length(phy);
n_ip = 2;
n_out = 2;
samp_T = 4;
Ns = 600;
N1 = 200;
N2 = 2 * N1 ;  


% for non linear simulation 
Xk = zeros(2,Ns); %state matrix
Tk = zeros(n_out,Ns);

% for estimation 
Tk_pred = zeros(n_out,Ns);
xk_pred = zeros(n_states,Ns);
yk_pred = zeros(n_states,Ns);
Uk_seq = zeros(n_ip,Ns);
error_k = zeros(n_out,Ns);
ek_f = zeros(n_out,Ns);
ref = zeros(n_out,Ns);
time = zeros(Ns,1);

%Weighing on Matrices 
Wx = C_mat'*C_mat;
Wu = eye(2,2);

%Bound on inputs
U_L = [0,0]';
U_H = [100,100]';
u_low = U_L - Us;
u_high = U_H - Us;


% Initialize controller tuning parameters
phy_e = diag( [ 0.98 0.98 ]) ;                 % Innovation filter parameter 
alfa_r = 0.9 ;                                       
phy_r = alfa_r * eye(n_out) ;                  % Setpoint filter parameter

%for simulations 
Ys = Ys';
Xk(:,1)= [29.1789 29.6676]' + [273.15 273.15]';
Tk(:,1) = Ys;
Uk_seq(:,1) = Us;
rk  = zeros(n_out,1);
rk = [57 57]' - Ys;
for k = 1:Ns
    time(k) = (k-1)*samp_T;
    %getting the measurements 
    % Specify inputs at instant k 
    if k == 1
        TCL.Uk(1)=Uk_seq(1,1);  
        TCL.Uk(2)=Uk_seq(2,1);
    else 
        TCL.Uk(1)=Uk_seq(1,k-1);  
        TCL.Uk(2)=Uk_seq(2,k-1);
    end
    % Process Simulation from sampling instant k to (k+1)
    [t,Xt]=ode45('TCL_Dynamics',[0 TCL.Samp_T], Xk(:,k));  
    Xk(:,k+1)=Xt(end,:)';
    Tk(:,k) = TCL.C_mat*Xk(:,k) - [273.15 273.15]' ;
    yk = Tk(:,k) - Ys ;  % Find perturbation measurement 
    if ( k <= (48/samp_T) )
        uk = [ 0 ; 0 ] ;   % Manual mode operation
    else
        if ( k < N1 )
            setpt = [ 57 57 ]' - Ys ;
        elseif (( k>= N1) && ( k < N2 ) )
            setpt = [65 50]' - Ys ;
        else
            setpt = [50 65]' - Ys ;
        end 
        % Generate filtered setpoint
        rk = phy_r * rk + (eye(n_out) - phy_r) * setpt ;
        
        % Compute target input and target states
        usk = inv_Ku_mat  * ( rk - Ke_mat * ek_f(:,k) ) ;
        xsk = inv_mat * ( gama * usk +  L_inf * ek_f(:,k) ) ;

        %%%%%%%%% MPC QP SOLVER %%%%%%%%%%%
        Wx = C_mat'*C_mat;
        Wu = eye(2,2);
        [Sx_mat, Su_mat, Se_mat, SxI_mat, SuI_mat, Wx_QP, Wu_QP ] = Generate_QP_Matrices( phy, gama, L_inf, Wx, Wu, N_pred, N_con);

        % Getting H and F matrices for QUADPROG 
        H = 2*(Su_mat'*Wx_QP*Su_mat+Wu_QP);
        disp((H))
        etak = Sx_mat*xk_pred(:,k) + Se_mat*ek_f(:,k) - SxI_mat*xsk;
        F = 2*(etak'*Wx_QP*Su_mat - (SuI_mat*usk)'*Wu_QP)';
        disp(length(F))

        % Getting A and B matrices for QUADPROG 
        A = [eye(n_ip*N_con,n_ip*N_con);-eye(n_ip*N_con,n_ip*N_con)];
        B = [SuI_mat*u_high;-SuI_mat*u_low];

        [U_f,optimum_cost] = quadprog(H,F,A,B);
        uk = U_f(1:2);
    end 
    Uk_seq(:,k) = uk + Us   ;
    ref(:,k) = rk + Ys ; 
          
    %h1(Uk(1));  % Send  manipulated inputs to heaters 
    %h2(Uk(2));
    
    % State Estimator computations from instant (k) to (k+1) 
     ek =  yk - C_mat *  xk_pred(:,k) ;    % Compute Innovation 
     xk_pred(:,k+1) = phy * xk_pred(:,k) + gama * uk + L_inf* ek_f(:,k)  ; 
     
    % Filtered innovation for innovation bias implementation
     ek_f(:,k+1) = phy_e * ek_f(:,k) + (eye(n_out) - phy_e) *  ek ;
     error_k(:,k) = ek;
end
figure(1), subplot(2,1,1) 
plot(time,Tk(1,:),'r-', time,ref(1,:),'k.-','MarkerSize',10), grid
ylabel('Temperature (degC)')
legend('Temperature 1','Setpoint 1')
subplot(2,1,2)
plot(time,Tk(2,:),'b-', time,ref(2,:),'k.-','MarkerSize',10), grid
legend('Temperature 2','Setpoint 2')
ylabel('Temperature (degC)')
xlabel('Time (sec)')

figure(2), subplot(2,1,1),
stairs(time,Uk_seq(1,:),'r-', 'LineWidth',2), grid
ylabel('Heater 1 (%)')
subplot(2,1,2), stairs(time,Uk_seq(2,:),'b-','LineWidth',2), grid
ylabel('Heater 2 (%)')
xlabel('Time (sec)')

figure(3), subplot(2,1,1) 
plot(time,error_k(1,:),'r-','MarkerSize',10), grid
title('Innovation Sequence') 
ylabel('ek_1(k)')
subplot(2,1,2) 
plot(time,error_k(2,:),'b-','MarkerSize',10), grid 
ylabel('ek_2(k)')
xlabel('Time (sec)')


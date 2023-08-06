% This program implements LQOC controller on TCL system using state space
% realization of an ARMAX 2th model identified from input-output data obtained on 
% 14 April 2022 
clear all 
clc 
%load Model_matrices.mat
load PRBS_Data_14_April_2022.mat 
load TCL_SS_Model_AMX_24July2020.mat
load TCL_LQOC_AMX_24July2020.mat
load TCL_Model_Parameters.mat%Experimental data obtained 
global TCL
Ys = mean(Yk_data(230:250,:));



sampl_T = 4;
samp_T = 4 ;
N_samp = 2400 /samp_T;
N1 = 800/samp_T ;
N2 = 2 * N1 ;  

n_st = length( phy ) ;  % Dimension of state variables 
n_ip = 2 ;                     % No. of inputs 
n_op = 2 ;                    % No. of outputs 

% Initialize controller tuning parameters
phy_e = diag( [ 0.98 0.98 ]) ;                 % Innovation filter parameter 
alfa_r = 0.8;                                       
phy_r = alfa_r * eye(n_op) ;                % Setpoint filter parameter 

xkpred = zeros(n_st,N_samp) ;                      % Initialize observer related variables 
uk = zeros(n_ip,1) ; 
ek_f = zeros( n_op, 1 ) ; 
Tk = Ys ; 
Uk = Us ; 
rk = zeros(n_op, 1) ;


t1s = zeros(N_samp,1);
t2s = zeros(N_samp,1);
h1s = zeros(N_samp,1);
h2s = zeros(N_samp,1);
e1s = zeros(N_samp,1);
e2s = zeros(N_samp,1) ;
R1s = zeros(N_samp,1);
R2s = zeros(N_samp,1);

Ys = Ys';    % Needed because Ys was saved as 1 x 2 instead of 2 x 1 vector 
Tk = zeros( n_op, 1) ;
rk = [ 57 57 ]' - Ys ;
time = zeros(N_samp,1) ;

Xk = zeros(2,N_samp);
Xk(:,1)= Yk_data(1,:)' + [273.15 273.15]';
for k = 1:N_samp
    tic;
    k
    time(k) = (k-1)*samp_T ;
    % Specify inputs at instant k 
    if k == 1  
        TCL.Uk(1)=h1s(1);  
        TCL.Uk(2)=h2s(1);  
    else 
        TCL.Uk(1)=h1s(k-1);  
        TCL.Uk(2)=h2s(k-1);
    end 
    % Process Simulation from sampling instant k to (k+1)
    [t,Xt]=ode45('TCL_Dynamics',[0 TCL.Samp_T], Xk(:,k));  
    Xk(:,k+1)=Xt(end,:)';
    Tk = TCL.C_mat*Xk(:,k) - [273.15 273.15]' ;
    yk = Tk - Ys ;  % Find perturbation measurement  

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
        rk = phy_r * rk + (eye(n_op) - phy_r) * setpt ;
        
        % Compute target input and target states
        usk = inv_Ku_mat  * ( rk - Ke_mat * ek_f ) ;
        xsk = inv_mat * ( gama * usk +  L_inf * ek_f ) ;
        
        % Innovation Bias controller implementation
        uk = usk  - G_inf * (xkpred(:,k) - xsk)  ;
        
        if ( uk(1) <= uk_L(1) )   % Impose input constraints on uk(1)
            uk(1) = uk_L(1) ;
        elseif ( uk(1) >= uk_H(1) )
            uk(1) = uk_H(1) ;
        end
        if ( uk(2) <= uk_L(2) ) % Inpose output constraints on uk(2)
            uk(2) = uk_L(2) ;
        elseif ( uk(2) >= uk_H(2) )
            uk(2) = uk_H(2) ;
        end
    end 
    
    Uk = uk + Us   ;
    Rk = rk + Ys ; 
          
    %h1(Uk(1));  % Send  manipulated inputs to heaters 
    %h2(Uk(2));
    
    % State Estimator computations from instant (k) to (k+1) 
     ek =  yk - C_mat *  xkpred(:,k) ;    % Compute Innovation 
     xkpred(:,k+1) = phy * xkpred(:,k) + gama * uk + L_inf* ek  ; 
     
    % Filtered innovation for innovation bias implementation
     ek_f = phy_e * ek_f + (eye(n_op) - phy_e) *  ek ;
    
    % plot heater and temperature data
    h1s(k) = Uk(1); h2s(k) = Uk(2) ;
    t1s(k) = Tk(1) ; t2s(k) = Tk(2) ;
    R1s(k) = Rk(1) ; R2s(k) = Rk(2)  ;
    e1s(k) = ek(1); e2s(k) = ek(2);
    
    
    t = toc;       
    %samp_T-t ;   % Time needed to implement control law and update figures 
    %pause(max(0.01, samp_T-t))     % Wait till next sampling interval 
end

%disp('Turn off heaters')
%h1(0);
%h2(0);

%disp('LQOC Test Complete')

% Plot results 

figure(1), subplot(2,1,1) 
plot(time,t1s,'r-', time,R1s,'k.-','MarkerSize',10), grid
axis([0 N_samp*samp_T 30 70])
ylabel('Temperature (degC)')
legend('Temperature 1','Setpoint 1')
subplot(2,1,2)
plot(time,t2s,'b-', time,R2s,'k.-','MarkerSize',10), grid
axis([0 N_samp*samp_T 30 70]) 
legend('Temperature 2','Setpoint 2')
ylabel('Temperature (degC)')
xlabel('Time (sec)')

figure(3), subplot(2,1,1),
stairs(time,h1s,'r-', 'LineWidth',2), grid
axis([0 N_samp*samp_T 0 100]) 
ylabel('Heater 1 (%)')
subplot(2,1,2), stairs(time,h2s,'b-','LineWidth',2), grid
axis([0 N_samp*samp_T 0 100]) 
ylabel('Heater 2 (%)')
xlabel('Time (sec)')

figure(4), subplot(2,1,1) 
plot(time,e1s,'r-','MarkerSize',10), grid
axis([0 N_samp*samp_T -5 5]) 
title('Innovation Sequence') 
ylabel('ek_1(k)')
subplot(2,1,2) 
plot(time,e2s,'b-','MarkerSize',10), grid
axis([0 N_samp*samp_T -5 5]) 
ylabel('ek_2(k)')
xlabel('Time (sec)')
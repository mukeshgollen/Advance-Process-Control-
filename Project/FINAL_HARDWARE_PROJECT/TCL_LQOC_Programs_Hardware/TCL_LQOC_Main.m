% This program implements LQOC controller on TCL system using state space
% realization of an ARMAX model identified from input-output data obtained on 
% 24 July 2022 

close all; clear all; clc

SetGraphics

tclab_N;       % Execute tclab_N.m for initialization of TCL system

disp('Load Model and Controller Files')
load TCL_LQOC_AMX_24July2020.mat   % This file loads (phy, gama, C,L_inf) and (Ys, Us)
load TCL_SS_Model_AMX_24July2020.mat     % This file loads LQOC related matrices 

samp_T = 4 ;
N_samp = 2400 /samp_T;
N1 = 800/samp_T ;
N2 = 2 * N1 ;  

n_st = length( phy ) ;  % Dimension of state variables 
n_ip = 2 ;                     % No. of inputs 
n_op = 2 ;                    % No. of outputs 

% Initialize controller tuning parameters
phy_e = diag( [ 0.9 0.98 ]) ;                 % Innovation filter parameter 
alfa_r = 0.85 ;                                       
phy_r = alfa_r * eye(n_op) ;                % Setpoint filter parameter 

xkpred = zeros(n_st,1) ;                      % Initialize observer related variables 
uk = zeros(n_ip,1) ; 
ek_f = zeros( n_op, 1 ) ; 
Tk = Ys ; 
Uk = Us ; 
rk = zeros(n_op, 1) ;

figure(1)
t1s = zeros(N_samp,1);
t2s = zeros(N_samp,1);
h1s = zeros(N_samp,1);
h2s = zeros(N_samp,1);
e1s = zeros(N_samp,1);
e2s = zeros(N_samp,1) ;
R1s = zeros(N_samp,1);
R2s = zeros(N_samp,1);

Ys = Ys' ;    % Needed because Ys was saved as 1 x 2 instead of 2 x 1 vector 
Tk = zeros( n_op, 1) ;
rk = [ 57 57 ]' - Ys ;
time = zeros(N_samp,1) ; 

for k = 1:N_samp
    tic;
    k
    time(k) = (k-1)*samp_T ;
    % read temperatures
    Tk(1) = T1C();
    Tk(2) = T2C();
    
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
        uk = usk  - G_inf * (xkpred - xsk)  ;
        
        if ( uk(1) <= uk_L(1) )   % Inpose input constraints on uk(1)
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
          
    h1(Uk(1));  % Send  manipulated inputs to heaters 
    h2(Uk(2));
    
    % State Estimator computations from instant (k) to (k+1) 
     ek =  yk - C_mat *  xkpred     % Compute Innovation 
     xkpred = phy * xkpred + gama * uk + L_inf* ek  ; 
     
    % Filtered innovation for innovation bias implementation
     ek_f = phy_e * ek_f + (eye(n_op) - phy_e) *  ek
    
    % plot heater and temperature data
    h1s(k) = Uk(1); h2s(k) = Uk(2) ;
    t1s(k) = Tk(1) ; t2s(k) = Tk(2) ;
    R1s(k) = Rk(1) ; R2s(k) = Rk(2)  ;
    e1s(k) = ek(1); e2s(k) = ek(2);
    
    if ( rem(k,3)==0)  % Update plot after every 15 seconds 
        clf
        subplot(2,1,1)
        plot(time(1:k),t1s(1:k),'r-', time(1:k),R1s(1:k),'k.-', 'LineWidth',2);
        hold on
        plot(time(1:k),t2s(1:k),'b-', time(1:k),R2s(1:k),'k.-', 'LineWidth',2);
        ylabel('Temperature (degC)')
        legend('Temperature 1','Temperature 2','Location','NorthWest')
        subplot(2,1,2)
        plot(time(1:k),h1s(1:k),'r-', 'LineWidth',2);
        hold on
        plot(time(1:k),h2s(1:k),'b-','LineWidth',2);
        ylabel('Heater (0-5.5 V)')
        xlabel('Time (sec)')
        legend('Heater 1','Heater 2','Location','NorthWest')
        drawnow;
    end
    t = toc;       
    samp_T-t    % Time needed to implement control law and update figures 
    pause(max(0.01, samp_T-t))     % Wait till next sampling interval 
end

disp('Turn off heaters')
h1(0);
h2(0);

disp('LQOC Test Complete')

% Plot results 

figure(2), subplot(2,1,1) 
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
      
save Results_LQOC_ARMX_4 time h1s h2s t1s t2s R1s R2s e1s e2s phy_e 
  
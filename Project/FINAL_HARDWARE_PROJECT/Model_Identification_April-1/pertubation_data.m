% The data collected from the TCL is used in order to have a pertubation
% data around the steady state for system identification and validation
clear all
load PRBS_Data_1_April_2022.mat %Experimental data 
Ys = mean(Yk_data(230:250,:));

% Generating Perturbation data for Model Identification 
a1 = 230;
a2 = 1000; 
Y_ident = Yk_data(a1:a2,:) -Ys;
U_ident = Uk_data(a1:a2,:) - Us';

TCL_identification =  iddata(Y_ident,U_ident,samp_T);
figure(1)
plot(TCL_identification),title('Identification-data'); 

% Generating Perturbation data for Model Validation 
a3 = 700;
Y_validate = Yk_data(a3:end,:) - Ys;
U_validate = Uk_data(a3:end,:) - Us';

TCL_validation = iddata(Y_validate,U_validate,samp_T);
figure(2)
plot(TCL_validation),title('Validation-data');


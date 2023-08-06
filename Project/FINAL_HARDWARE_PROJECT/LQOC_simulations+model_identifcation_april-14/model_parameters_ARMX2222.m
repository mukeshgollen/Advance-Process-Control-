%let's calculate the optimal gain matrix  for the system by solving the riccati
%equation 

W_u = [1,0;0,1];
%W_x = [5,0,0,0,0,0,0,0;0,10,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0];
W_y = [5,0;0,10];

[ phy_m, gama_m, C_mat_m, D_mat_m, L_mat_m, x0_m] = idssdata(ss1);
W_x = C_mat_m'*W_y*C_mat_m;
[G_inf,S_inf,E_v] = dlqr(phy_m,gama_m,W_x,W_u);
 
inv_mat = inv(eye(2,2) - phy_m);
m = inv_mat\gama_m;
Ku_mat = C_mat_m*m;
n = inv_mat\L_mat_m;
Ke_mat = C_mat_m*n + eye(2,2);
inv_Ku_mat = inv(Ku_mat);
uk_L = [-30;-40];
uk_H = [70;60];
gamma = gama_m;
phy = phy_m;
C_mat = C_mat_m;
L_inf = L_mat_m;
sampl_T = 4;
Ys = mean(Yk_data(230:250,:));
save modelparamters_LQOC_state_space_new W_x W_u G_inf inv_mat Ke_mat Ku_mat inv_Ku_mat uk_L uk_H
save SS_Model_state_space_new  C_mat L_inf gamma phy sampl_T Ys Us ss1
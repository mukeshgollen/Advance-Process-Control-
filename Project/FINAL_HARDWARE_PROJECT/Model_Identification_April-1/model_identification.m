% Using the System Identification Toolbox we finally have our 4th order ARX
% model. Now we will obtain the State space realization of the ARX Model 

[ phy_model, gama_model, C_mat_model, D_mat_model, L_mat_model, x0_model] = idssdata( amx4444);
Modelfromdata.phy_m = phy_model;
Modelfromdata.gama_m = gama_model;
Modelfromdata.C_mat_m = C_mat_model;
Modelfromdata.D_mat_m = D_mat_model;
Modelfromdata.L_mat_m = L_mat_model;
save Model_matrices Modelfromdata
save ARMX4444 amx4444
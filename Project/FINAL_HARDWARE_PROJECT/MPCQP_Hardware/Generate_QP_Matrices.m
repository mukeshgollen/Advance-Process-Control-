% This function generates matrices for reformulating MPC as a QP 
% Returns a structure mpc that contains these matrices

% This function generates block matrices. While genarting the block
% matrices, the last row is compatible with the terminal state

function [Sx_mat, Su_mat, Se_mat, SxI_mat, SuI_mat, Wx_QP, Wu_QP ] = Generate_QP_Matrices( phy, gama, Lp_mat, Wx, Wu, N_pred, N_con)
  
% Innovation bias approach 

[n_st, n_ip] = size(gama) ;
n_op = n_st ; 
C_mat = eye(n_st) ; 
Gama_beta = Lp_mat ;
C_eta = eye(n_op) ;

% Intermediate matrices for calculations 
Sx_mat = [] ;     
Su_mat = [] ;
SxI_mat = [] ;     
SuI_mat = [] ;
Sbeta_mat = [];
Seta_mat= [];
Wx_QP = [] ;
Wu_QP = [] ;

mat = zeros(n_st, n_ip*N_pred);
c2 = 0 ;
temp_mat = eye(size(phy)) ;

for k = 1 : N_pred
    
    Wx_QP = blkdiag( Wx_QP, Wx ) ;
    SxI_mat = vertcat( SxI_mat, eye(n_st) ) ;  
    
    % Generation of Seta_mat
    Seta_mat = [Seta_mat; C_eta] ;

    % Generation of Sbeta_mat
    Sbeta_mat = [Sbeta_mat ; C_mat*temp_mat*Gama_beta];
    if ( k < N_pred )
        temp_mat = phy*temp_mat + eye(size(phy)) ;
    end
    
    mat = phy*mat;
    
    % Generation of Su_mat   
    c1 = c2 ; c2 = c1 + n_ip ;
    mat(:,c1+1:c2) = gama;    
    Su_mat = [ Su_mat ; C_mat*mat ] ;
    
    % Generation of Sx_mat
    Sx_mat = [ Sx_mat ; C_mat*phy^k ] ;
end

Se_mat = Sbeta_mat ; 

%  Generation of transformation mattrix for restricting the input move to control horizon
psi_pq_mat = [];
for i = 1:N_con-1
    psi_pq_mat = blkdiag( psi_pq_mat, eye(n_ip) );
end
I_mi = [];
for j = 1:(N_pred-N_con + 1)
    I_mi = vertcat( I_mi, eye(n_ip) ) ;
end
psi_pq_mat = blkdiag( psi_pq_mat, I_mi );

% Generation of Su_mat taking into account control horizon
Su_mat = Su_mat * psi_pq_mat;

  for k = 1 : N_con
        Wu_QP = blkdiag( Wu_QP, Wu ) ;
        SuI_mat = vertcat( SuI_mat, eye(n_ip) ) ;  
 end

return

end 


function fofX = System1_Dynamics(X,M,N,alpha)


fofX = zeros(3,1); 

z = alpha(1) * exp(-alpha(2)/X(2));

fofX(1) = (M(1)/100)*(N-X(1)) - 2*z*X(1)^2;
fofX(2) = (M(1)/100)* (275-X(2)) + alpha(3)*z*X(1)^2 - alpha(4)*(X(2)-X(3)) ;
fofX(3) = (M(2)/10)*(250 -X(3)) + alpha(5)*(X(2)-X(3));
end 


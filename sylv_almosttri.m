% sylv_almosttri.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
% Solves the Sylvester equation
% $\mathbf A\mathbf X+\mathbf X\mathbf B=\mathbf C$,
% with $\mathbf A$ being almost lower triangular, except for $a_{12}\neq0$.
function X=sylv_almosttri(A,B,C)
nA=size(A,1); % Get the order of $\mathbf A$
nB=size(B,1); % Get the order of $\mathbf B$
X=zeros(nA,nB); % Preallocate $\mathbf X$
% Obtain the Schur decomposition of $\mathbf A(1:2,1:2)$
a11=A(1,1); a12=A(1,2); a21=A(2,1); a22=A(2,2);
lambda1=((a11+a22)+sqrt((a11-a22)^2+4*a21*a12))/2; % Eigenvalue $\lambda_1$
Q=[a12,conj(a11-lambda1);lambda1-a11,conj(a12)]...
    /sqrt(abs(a12)^2+abs(lambda1-a11)^2); % Unitary matrix $\mathbf Q$
T=Q'*A(1:2,1:2)*Q; % Upper triangular matrix $\mathbf T$
tildeC=Q'*C(1:2,:);
tildex2=tildeC(2,:)/(T(2,2)*eye(nB)+B);
tildex1=(tildeC(1,:)-T(1,2)*tildex2)/(T(1,1)*eye(nB)+B);
X(1:2,:)=Q*[tildex1;tildex2]; % Obtain the first rows of $\mathbf X$
for i=3:nA % Obtain the remaining rows of $\mathbf X$.
    X(i,:)=(C(i,:)-A(i,1:i-1)*X(1:i-1,:))/(A(i,i)*eye(nB)+B);
end
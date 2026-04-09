% TestSylvester.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
% Compares the solvers sylv_almosttri, lyap, and sylvester in terms of
% elapsed time and accuracy, when solving the Sylvester equation
% $\mathbf A\mathbf X+\mathbf X\mathbf B=\mathbf C$,
% with $\mathbf A$ being almost lower triangular, except for $a_{12}\neq0$.
clear
rng(1); % Fix the seed to obtain reproducible matrices
nA=5000; % Order of $\mathbf A$
nB=20; % Order of $\mathbf B$
A=randn(nA)+1i*randn(nA); % Generate $\mathbf A$
A=tril(A); % $\mathbf A$ is almost lower-triangular
A(1,2)=randn+1i*randn; % Set $a_{12}\neq0$.
B=randn(nB)+1i*randn(nB); % Generate $\mathbf B$
% Enforce that $\mathbf A$ and $\mathbf B$ are strictly diagonally dominant
A(1:nA+1:end)=0; % Set the diagonal to zero
A=A+diag(sum(abs(A),2)+1i*randn(nA,1));
B(1:nB+1:end)=0; % Set the diagonal to zero
B=B+diag(sum(abs(B),2)+1i*randn(nB,1));
X=randn(nA,nB)+1i*randn(nA,nB); % Generate $\mathbf X$
C=A*X+X*B; % Compute $\mathbf C$
tic, X1=sylv_almosttri(A,B,C); toc % Elapsed time of sylv_almosttri
tic, X2=lyap(A,B,-C); toc % Elapsed time of lyap
tic, X3=sylvester(A,B,C); toc % Elapsed time of sylvester
disp(norm(X(:)-X1(:),inf)) % Maximum error of sylv_almosttri
disp(norm(X(:)-X2(:),inf)) % Maximum error of lyap
disp(norm(X(:)-X3(:),inf)) % Maximum error of sylvester
% PDE1.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
clear
tic
a=0.17; % $\alpha$
Nx=16; % $N_x$
Nt=2700; % $N_t$
tf=1.2; % $t_f$
t=tf*(0:Nt)'/Nt; % $t$
Dta=GenerateDa(Nt,tf,a); % $\mathbf D_t^\alpha$
b=1.4; % Scale factor $b$
[x,DD]=herdif(Nx,2,b);
x=x.'; % $x$
Dx=DD(:,:,1).'; % $\mathbf D_x$
Dx2=DD(:,:,2).'; % $\mathbf D_x^2$
[T,X]=ndgrid(t,x); % Two-dimensional grid
u=@(t,x)exp(2*t-x.^2); % $u(t, x)$
u0=@(x)exp(-x.^2); % $u_0(x)$
Et=[zeros(1,Nt);eye(Nt)]; % $\mathbf E_t$
Ft=[u0(x);zeros(Nt,Nx)]; % $\mathbf F_t$
A2=diag(2*x); % $\mathbf A_2$
A3=2*eye(Nx); % $\mathbf A_3$
A4=(2^a/gamma(1-a))*(gamma(1-a)-igamma(1-a,2*t)).*exp(2*T-X.^2); % $\mathbf A_4$
Bx=Dx2+Dx*A2+A3; % $\mathbf B_x$
A=Et'*Dta*Et; % $\mathbf A$
B=-Bx; % $\mathbf B$
C=-Et'*Dta*Ft+Et'*A4; % $\mathbf C$
Uinner=SylvesterAAlmostTriangular(A,B,C); % $\mathbf U_{inner}$
Unum=Et*Uinner+Ft; % $\mathbf U_{num}$
disp(max(max(abs(Unum-u(T,X))))) % Error
toc % Elapsed time

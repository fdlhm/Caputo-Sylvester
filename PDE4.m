% PDE4.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
clear
tic
a=0.15; % $\alpha$
Nt=3000; % $N_t$
Nx=10; % $N_x$
tf=1; % $t_f$
t=tf*(0:Nt)'/Nt; % $t$
Dta=GenerateDa(Nt,tf,a); % $\mathbf D_t^\alpha$
[Dx,x]=cheb(Nx);
x=(x+1)'/2; % $x$
Dx=2*Dx'; % $\mathbf D_x$
[T,X]=ndgrid(t,x); % Two-dimensional grid
u=@(t,x)exp(x).*sqrt(t); % $u(t, x)$
ua=@(t)sqrt(t); % $u_a(t)$
ub=@(t)exp(1)*sqrt(t); % $u_b(t)$
Et=[zeros(1,Nt);eye(Nt)]; % $\mathbf E_t$
Ex=[zeros(Nx-1,1),eye(Nx-1),zeros(Nx-1,1)]; % $\mathbf E_x$
Fx=[ub(t),zeros(Nt+1,Nx-1),ua(t)]; % $\mathbf F_x$
A4=((sqrt(pi)/2/gamma(3/2-a))*T.^(1/2-a)-sqrt(T)).*exp(X); % $\mathbf A_4$
Bx=Dx^2-Dx+eye(Nx+1); % $\mathbf B_x$
A=Et'*Dta*Et; % $\mathbf A$
B=-Ex*Bx*Ex.'; % $\mathbf B$
C=(-Et'*Dta*Et*Fx(2:end,:)+Fx(2:end,:)*Bx+Et'*A4)*Ex.'; % $\mathbf C$
Uinner=SylvesterAAlmostTriangular(A,B,C); % $\mathbf U_{inner}$
Unum=Et*Uinner*Ex+Et*Fx(2:end,:); % $\mathbf U_{num}$
disp(max(max(abs(Unum-u(T,X))))) % Error
toc % Elapsed time
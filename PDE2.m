% PDE2.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
clear
tic
a=0.17; % $\alpha$
Nt=2700; % $N_t$
Nx=15; % $N_x$
xa=-1.1; % $x_a$
xb=1.3; % $x_b$
tf=1.2; % $t_f$
t=tf*(0:Nt)'/Nt; % $t$
Dta=GenerateDa(Nt,tf,a); % $\mathbf D_t^\alpha$
[Dx,x]=cheb(Nx);
x=xa+(xb-xa)*(x+1)'/2; % $x$
Dx=(2/(xb-xa))*Dx'; % $\mathbf D_x$
[T,X]=ndgrid(t,x); % Two-dimensional grid
u=@(t,x)exp(2*t+1.5*x); % $u(t, x)$
u0=@(x)exp(1.5*x); % $u_0(x)$
ca=1;da=2; % $c_a$ and $d_a$
ua=@(t)4*exp(2*t-1.65); % $u_a(t)$
cb=3;db=4; % $c_b$ and $d_b$
ub=@(t)9*exp(2*t+1.95); % $u_b(t)$
Et=[zeros(1,Nt);eye(Nt)]; % $\mathbf E_t$
Ft=[u0(x);zeros(Nt,Nx+1)]; % $\mathbf F_t$
c11=da*Dx(1,Nx+1);c12=ca+da*Dx(Nx+1,Nx+1); % $c_{11}$ and $c_{12}$
c21=cb+db*Dx(1,1);c22=db*Dx(Nx+1,1); % $c_{21}$ and $c_{22}$
Ex=[(-da*c22*Dx(2:Nx,Nx+1)+db*c12*Dx(2:Nx,1))/(c11*c22-c12*c21),eye(Nx-1),...
    (da*c21*Dx(2:Nx,Nx+1)-db*c11*Dx(2:Nx,1))/(c11*c22-c12*c21)]; % $\mathbf E_x$
Fx=[(c22*ua(t)-c12*ub(t))./(c11*c22-c12*c21),zeros(Nt+1,Nx-1),...
    (-c21*ua(t)+c11*ub(t))./(c11*c22-c12*c21)]; % $\mathbf F_x$
A1=diag((2^a/2.25)*(1+x.^2)); % $\mathbf A_1$
A2=diag((2^a/1.5)*x.^2); % $\mathbf A_2$
A3=diag(-2^(a+1)*x.^2); % $\mathbf A_3$
A4=-(2^a*igamma(1-a,2*t)/gamma(1-a)).*exp(2*T+1.5*X); % $\mathbf A_4$
Bx=Dx^2*A1+Dx*A2+A3; % $\mathbf B_x$
A=Et'*Dta*Et; % $\mathbf A$
B=-Ex*Bx/Ex; % $\mathbf B$
C=(-Et'*Dta*Et*Fx(2:end,:)-Et'*Dta*Ft+Fx(2:end,:)*Bx+Et'*A4)/Ex; % $\mathbf C$
Uinner=sylv_almosttri(A,B,C); % $\mathbf U_{inner}$
Unum=Et*Uinner*Ex+Et*Fx(2:end,:)+Ft; % $\mathbf U_{num}$
disp(max(max(abs(Unum-u(T,X))))) % Error
toc % Elapsed time
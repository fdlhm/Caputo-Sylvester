% CaputoDerivative.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
clear
tic
a=0.17; % $\alpha$
tf=1.2; % $t_f$
N=100; % $N$
tj=tf*(0:N)'/N; % $t_j$
h=tf/N; % $h$
f=exp(2*tj); % $f(t_j)$
Dtaf=(2^a/gamma(1-a))*(exp(2*tj).*(gamma(1-a)-igamma(1-a,2*tj))); % $D_t^\alpha f$
toc,tic % Elapsed time needed to generate $D^\alpha f$
Dtafnum=zeros(N+1,1); % $D_{t,num}^\alpha f$
Dtafnumstar=zeros(N+1,1); % $D_{t,num}^{\alpha,*}f$
aux1=(0:N).^(1-a); % Precompute $j^{1-\alpha}$
aux2=(0:N).^(2-a); % Precompute of $j^{2-\alpha}$
for j=2:N
    Dtafnum(j+1)=Dtafnum(j+1)...
        +(aux2(j:-1:2)-aux2(j-1:-1:1))*(f(3:j+1)-2*f(2:j)+f(1:j-1))/(2-a)...
        +.5*aux1(j:-1:2)*(f(3:j+1)-f(1:j-1))...
        -.5*aux1(j-1:-1:1)*(3*f(3:j+1)-4*f(2:j)+f(1:j-1));
end
Dtafnumstar(2:N+1)=Dtafnum(2:N+1)+(f(2)-f(1))*(aux1(2:N+1)-aux1(1:N))';
Dtafnum(2:N+1)=Dtafnum(2:N+1)+(f(3)-2*f(2)+f(1))/(2-a)*(aux2(2:N+1)-aux2(1:N))'...
    -.5*(f(3)-4*f(2)+3*f(1))*aux1(2:N+1)'-.5*(f(3)-f(1))*aux1(1:N)';
Dtafnum=(h^(-a)/gamma(2-a))*Dtafnum; % $D_{t,num}^\alpha f$
Dtafnumstar=(h^(-a)/gamma(2-a))*Dtafnumstar; % $D_{t,num}^{\alpha,*}f$
toc % Elapsed time needed to generate $D_{t,num}^\alpha f$ and $D_{t,num}^{\alpha,*}f$
disp(norm(Dtaf-Dtafnum,inf))
disp(norm(Dtaf-Dtafnumstar,inf))
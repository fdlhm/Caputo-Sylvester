% FastConvolution.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
clear
tic
a=0.17; % $\alpha$
tf=1.2; % $t_f$
N=2^17; % $N$
tj=tf*(0:N)'/N; % $t_j$
h=tf/N; % $h$
f=exp(2*tj); % $f(t_j)$
Dtaf=(2^a/gamma(1-a))*(exp(2*tj).*(gamma(1-a)-igamma(1-a,2*tj))); % $D_t^\alpha f$
toc,tic % Elapsed time needed to generate $D^\alpha f$
P=2^ceil(log2(2*N-3)); % Length $P$ of $\tilde c_1$ and the other vectors
aux1=(0:N).^(1-a); % Precompute $j^{1-\alpha}$
aux2=(0:N).^(2-a); % Precompute $j^{2-\alpha}$
a1tilde=[(f(3:N+1)-2*f(2:N)+f(1:N-1))/(2-a);zeros(P-N+1,1)]; % $\tilde a_1$
a2tilde=[.5*(f(3:N+1)-f(1:N-1));zeros(P-N+1,1)]; % $\tilde a_2$
a3tilde=[-.5*(3*f(3:N+1)-4*f(2:N)+f(1:N-1));zeros(P-N+1,1)]; % $\tilde a_3$
b1tilde=[(aux2(2:N)-aux2(1:N-1))';zeros(P-N+1,1)]; % $\tilde b_1$
b2tilde=[aux1(2:N)';zeros(P-N+1,1)]; % $\tilde b_2$
b3tilde=[aux1(1:N-1)';zeros(P-N+1,1)]; % $\tilde b_3$
ctilde=ifft(fft(a1tilde).*fft(b1tilde)+fft(a2tilde).*fft(b2tilde)...
    +fft(a3tilde).*fft(b3tilde)); % $\tilde c_1$
Dtafnum=[0;(f(3)-2*f(2)+f(1))/(2-a)-.5*(f(3)-4*f(2)+3*f(1));
    ((f(3)-2*f(2)+f(1))/(2-a))*(aux2(3:N+1)-aux2(2:N))'...
        - .5*(f(3)-4*f(2)+3*f(1))*(2:N)'.^(1-a)...
        - .5*(f(3)-f(1))*aux1(2:N)'+ctilde(1:N-1)];
Dtafnum=(h^(-a)/gamma(2-a))*Dtafnum; % $D_{t,num}^\alpha f$
toc % Elapsed time needed to generate $D_{t,num}^\alpha f$
disp(norm(Dtaf-Dtafnum,inf))
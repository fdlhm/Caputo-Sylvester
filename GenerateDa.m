% GenerateDa.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via Sylvester equations,
% arXiv:2501.09180, (2025).
%
function Dta=GenerateDa(N,tf,a)
Dta=zeros(N+1); % Preallocate $\mathbf D_t^\alpha$
aux1=(0:N).^(1-a); % Precompute $j^{1-\alpha}$
aux2=(0:N).^(2-a); % Precompute $j^{2-\alpha}$
Dta(2:N+1,1:3)=[(aux2(2:N+1)-aux2(1:N))/(2-a)-1.5*aux1(2:N+1)+.5*aux1(1:N);
    -(2/(2-a))*(aux2(2:N+1)-aux2(1:N))+2*aux1(2:N+1);
    (aux2(2:N+1)-aux2(1:N))/(2-a)-.5*aux1(2:N+1)-.5*aux1(1:N)]';
for j=3:N+1
    Dta(j,1:j-2)=Dta(j,1:j-2)+(aux2(j-1:-1:2)-aux2(j-2:-1:1))/(2-a)...
        -.5*aux1(j-1:-1:2)-.5*aux1(j-2:-1:1);
    Dta(j,2:j-1)=Dta(j,2:j-1)-(2/(2-a))*(aux2(j-1:-1:2)-aux2(j-2:-1:1))...
        +2*aux1(j-2:-1:1);
    Dta(j,3:j)=Dta(j,3:j)+(aux2(j-1:-1:2)-aux2(j-2:-1:1))/(2-a)...
        +.5*aux1(j-1:-1:2)-1.5*aux1(j-2:-1:1);
end
Dta=((tf/N)^(-a)/gamma(2-a))*Dta;
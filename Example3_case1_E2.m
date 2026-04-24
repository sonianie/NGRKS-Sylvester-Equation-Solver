% A,B,F,G come from SC2024 -- Example 2
% A. Casulli and L. Robol, An efficient block rational Krylov solver for Sylvester equations with
% adaptive pole selection, SIAM J. Sci. Comput., 46 (2024), A798–A824.

clc;
clear;

% compute the approximate solution for the convection-diffusion partial differential equation (SC2024, Example 2)
n = 2^11; % Size of the problem
[A,B,U,V]=matrices_example2(n); V = -V;

% Set up solver parameters including tolerance and estimated eigenvalue ranges
As2  = norm(A,1);
As1  = As2/condest(A);
Bs2  = norm(B,1);
Bs1  = Bs2/condest(B);
tol  = 1e-8;

% Run NGRKS(1) algorithm and record performance matrics
tic;
[~,~,itG,RESG] = NGRKSn_sylv(A,-B,-U,V,As1,As2,Bs1,Bs2,1,tol);
cpu(1)=toc;
it(1)=itG;
residual(1)=RESG(end);
ratio(1) = cpu(1)/it(1);

% Run NGRKS(4) algorithm and record performance matrics
tic;
[~,~,itG4,RESG4] = NGRKSn_sylv(A,-B,-U,V,As1,As2,Bs1,Bs2,4,tol);
cpu(2)=toc;
it(2)=itG4;
residual(2)=RESG4(end);
ratio(2) = cpu(2)/it(2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,U,V]=matrices_example2(n)
h=1/(n+1);
vi=0.0083;
w1=@(x)(1+(1+x).^2/4);
w2=@(y) (1/2)*y;

t=linspace(0,1,n).';
Phi1=spdiags(w1(t), 0, n, n);
Psi2=spdiags(w2(t), 0, n, n);

[X,Y]=meshgrid(t);
F=1./((1+X+Y));
[U,S,V] = svd(F);
c=2;
while c>1
    if S(c,c)<1e-10
        S=S(1:c-1,1:c-1);
        U=U(:,1:c-1);
        V=V(:,1:c-1);
        c=0;
    end
    c=c+1;
end
U=U*sqrt(S);
V=V*sqrt(S);
T1=spdiags(ones(n,1) * [-1 2 -1],-1:1, n, n);
T1=(1/h^2)*T1;
T2=T1';
B1=spdiags(ones(n,1) * [-1 1], [-1 1], n, n);
B1=(1/(2*h))*B1;
B2=B1';
A=-vi*T1+Phi1*B1;
B=+vi*T2-B2*Psi2;
end

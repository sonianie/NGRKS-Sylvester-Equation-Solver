clc;
clear;


%marices A and B is obtained from the discretisation of the operator
%Lu = ∆u − f1(x, y)∂u/∂x + f2(x, y)∂u/∂y + g(x, y)
%on the unit square [0, 1]×[0, 1] with homogeneous Dirichlet boundary conditions

n0 = 200; m0 = 100;
afx_str = 'cos(x+y)'; afy_str = 'sin(x-y)'; ag_str = 'x*y';
[A, ~] = fdm_2d_matrix(n0,afx_str,afy_str,ag_str);
As2  = norm(A,1);
As1  = As2/condest(A);

bfx_str = 'x+2*y'; bfy_str = 'exp(y-x)'; bg_str = 'x^2-y^2';
[B, ~] = fdm_2d_matrix(m0,bfx_str,bfy_str,bg_str);
Bs2  = norm(B,1);
Bs1  = Bs2/condest(B);

% Set up solver parameters including tolerance and estimated eigenvalue ranges
n    = size(A,2);  m = size(B,2);
tol  = 1e-8;
PP   = [];
CPU1 = [];IT1 = [];RA1 = [];
CPU2 = [];IT2 = [];RA2 = [];


for p = [10,20,30]
    PP = [PP,p];
    F = rand(n,p); G = rand(m,p); 

    % Run NGRKS(c) with c = 1,4 and record performance matrics
    tic;
    [~,~,it1,RES1] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,1,tol);
    cpu1 = toc;
    ratio1 = cpu1/it1;
    tic;
    [~,~,it2,RES2] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,4,tol);
    cpu2 = toc;
    ratio2 = cpu2/it2;


    CPU1 = [CPU1,cpu1];IT1 = [IT1,it1];RA1 = [RA1,ratio1];
    CPU2 = [CPU2,cpu2];IT2 = [IT2,it2];RA2 = [RA2,ratio2];
end


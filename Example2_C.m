clc;
clear;

% Generate 2D finite difference matrices A and B based on given function strings
afx_str = 'cos(x+y)'; afy_str = 'sin(x-y)'; ag_str = 'x*y';
[A, ~] = fdm_2d_matrix(300,afx_str,afy_str,ag_str);
As2  = norm(A,1);
As1  = As2/condest(A);

bfx_str = 'x+2*y'; bfy_str = 'exp(y-x)'; bg_str = 'x^2-y^2';
[B, ~] = fdm_2d_matrix(250,bfx_str,bfy_str,bg_str);
Bs2  = norm(B,1);
Bs1  = Bs2/condest(B);

% Set up problem dimensions and solver parameters
n = size(A,2);  m = size(B,2); cpu_max=500;  tol  = 1e-8;

% Initialize arrays to store performance metrics for different cycle lengths
CPU1 = [];IT1 = [];RR1 = [];RANK1 = [];
CPU2 = [];IT2 = [];RR2 = [];RANK2 = [];
CPU3 = [];IT3 = [];RR3 = [];RANK3 = [];
CPU4 = [];IT4 = [];RR4 = [];RANK4 = [];
CPU5 = [];IT5 = [];RR5 = [];RANK5 = [];
CPU6 = [];IT6 = [];RR6 = [];RANK6 = [];
CPU7 = [];IT7 = [];RR7 = [];RANK7 = [];
CPU8 = [];IT8 = [];RR8 = [];RANK8 = [];
CPU9 = [];IT9 = [];RR9 = [];RANK9 = [];

% Run experiments for different cycle lengths (c = 1 to 9) with fixed p=20
for p = 20
    F = rand(n,p); G = rand(m,p); 
    
    % Run NGRKS algorithm with different cycle parameter c and measure performance
    tic;[Z1,~,it1,RES1] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,1,tol);cpu1=toc;
    rankZ1 = rank(Z1);Z1 =[];
    tic;[Z2,~,it2,RES2] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,2,tol);cpu2=toc;
    rankZ2 = rank(Z2);Z2 =[];
    tic;[Z3,~,it3,RES3] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,3,tol);cpu3=toc;
    rankZ3 = rank(Z3);Z3 =[];
    tic;[Z4,~,it4,RES4] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,4,tol);cpu4=toc;
    rankZ4 = rank(Z4);Z4 =[];
    tic;[Z5,~,it5,RES5] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,5,tol);cpu5=toc;
    rankZ5 = rank(Z5);Z5 =[];
    tic;[Z6,~,it6,RES6] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,6,tol);cpu6=toc;
    rankZ6 = rank(Z6);Z6 =[];
    tic;[Z7,~,it7,RES7] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,7,tol);cpu7=toc;
    rankZ7 = rank(Z7);Z7 =[];
    tic;[Z8,~,it8,RES8] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,8,tol);cpu8=toc;
    rankZ8 = rank(Z8);Z8 =[];
    tic;[Z9,~,it9,RES9] = NGRKSn_sylv(A,B,F,G,As1,As2,Bs1,Bs2,9,tol);cpu9=toc;
    rankZ9 = rank(Z9);Z9 =[];

     % Store performance metrics (CPU time, iterations, final residual, rank) for each c
    CPU1 = [CPU1,cpu1];IT1 = [IT1,it1];RR1 = [RR1,RES1(end)];RANK1 =[RANK1,rankZ1];
    CPU2 = [CPU2,cpu2];IT2 = [IT2,it2];RR2 = [RR2,RES2(end)];RANK2 =[RANK2,rankZ2];
    CPU3 = [CPU3,cpu3];IT3 = [IT3,it3];RR3 = [RR3,RES3(end)];RANK3 =[RANK3,rankZ3];
    CPU4 = [CPU4,cpu4];IT4 = [IT4,it4];RR4 = [RR4,RES4(end)];RANK4 =[RANK4,rankZ4];
    CPU5 = [CPU5,cpu5];IT5 = [IT5,it5];RR5 = [RR5,RES5(end)];RANK5 =[RANK5,rankZ5];
    CPU6 = [CPU6,cpu6];IT6 = [IT6,it6];RR6 = [RR6,RES6(end)];RANK6 =[RANK6,rankZ6];
    CPU7 = [CPU7,cpu7];IT7 = [IT7,it7];RR7 = [RR7,RES7(end)];RANK7 =[RANK7,rankZ7];
    CPU8 = [CPU8,cpu8];IT8 = [IT8,it8];RR8 = [RR8,RES8(end)];RANK8 =[RANK8,rankZ8];
    CPU9 = [CPU9,cpu9];IT9 = [IT9,it9];RR9 = [RR9,RES9(end)];RANK9 =[RANK9,rankZ9];
    
end

% Plot the convergence history for all c = 1 to 9
hold on
box on
plot((1:it1),log10(RES1),'r-', 'LineWidth', 2.5);
plot((1:it2),log10(kron(RES2,ones(2,1))),'g--', 'LineWidth', 2.5);
plot((1:it3),log10(kron(RES3,ones(3,1))),'b:', 'LineWidth', 2.5);
plot((1:it4),log10(kron(RES4,ones(4,1))),'m-.', 'LineWidth', 2.5);
plot((1:it5),log10(kron(RES5,ones(5,1))),'c-o', 'MarkerSize', 5, 'LineWidth', 2);
plot((1:it6),log10(kron(RES6,ones(6,1))),'k-s', 'MarkerSize', 5, 'LineWidth', 2);
plot((1:it7),log10(kron(RES7,ones(7,1))),'y-d', 'MarkerSize', 5, 'LineWidth',2);
plot((1:it8),log10(kron(RES8,ones(8,1))),'-p', 'Color', [0.5 0.5 0], 'MarkerSize', 5, 'LineWidth',2);
plot((1:it9),log10(kron(RES9,ones(9,1))),'--*', 'Color', [0.5 0 0.5], 'MarkerSize', 5, 'LineWidth', 2);
xlabel('Iter');
ylabel('log_{10}(RES) ');
legend('c=1','c=2','c=3','c=4','c=5','c=6','c=7','c=8','c=9');



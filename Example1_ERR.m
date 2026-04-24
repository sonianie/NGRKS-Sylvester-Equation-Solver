clc;
clear;

% Load and set up the test matrices A and B from Matrix Market files
% the University of Florida Sparse Matrix Collection can be downed at https://sparse.tamu.edu
A = -mmread('nos6.mtx');
B = -mmread('nos5.mtx');
% Set up dimensions and generate random right-hand side matrices F and G
n = size(A,2); m = size(B,2); p=2;
F = rand(n,p); G = rand(m,p); 

% Set up solver parameters including tolerance and estimated eigenvalue ranges
tol  = 1e-8;
As2  = norm(A,1);
As1  = As2/condest(A);
Bs2  = norm(B,1);
Bs1  = Bs2/condest(B);


% Call the error computation function for the NGRKS method
[iter,ER1,ER2] = err_NGRKSn_sylv(A,B,F,G,1,As1,As2,Bs1,Bs2,tol);

% Plot the computed errors and their theoretical bounds
figure(1)
hold on
box on
plot((1:length(ER1)),log10(ER1),'r-', 'LineWidth', 2.5);
plot((1:length(ER2)),log10(ER2),'b*-', 'LineWidth', 2.5);
xlabel('Iter');
ylabel('log_{10}(Err)');
legend('real valus of Err', 'upper bounds given by (4.9)');



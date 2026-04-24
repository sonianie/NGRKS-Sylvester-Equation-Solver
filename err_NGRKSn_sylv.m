function [iter,E1,E2] = err_NGRKSn_sylv(A,B,F,G,km,As1,As2,Bs1,Bs2,tol)
% Error norm Estimation of Novel Global Rational Krylov Subspace method for Sylvester matrix equations: AX+XB+F*G'=0
% Input:
%       A \in R^{n \times n},     B \in R^{l \times l}
%       F \in R^{n \times p},     G \in R^{l \times p}
%       As1(Bs1):    the smallest eigenvalues of A(B)
%       As2(Bs2):    the largest eigenvalues of A(B)
%       km:          cyclic constant
%       tol:         tolerance
% Output
%       iter:        the number of iterations of NGRKS(km)
%       E1:          the true error norm of NGRKS(km)
%       E2:          the estimated error norm of NGRKS(km)
% niesong1995@nuaa.edu.cn 2022.02.23
%%

% --- Initialization Section ---                
n       = size(A,2);       m      = size(B,2);    p    = size(F,2);
HA      = [];              HB     = [];          
SSA     = [];              SSB    = [];      
KKAA0   = [];              KKBB0  = [];
E1      = [];              E2     = [];
iter    = 0; 

% Calculate the exact solution of Sylvester matrix equations AX+XB+F*G'=0
XX = lyap(A,B,F*G');

% Calculate norms of F and G, initialize starting vectors V and W
nrmF   = norm(F,'fro');        
Avv    = F/nrmF;           V   = Avv;                    
EEn    = speye(n);  
nrmG   = norm(G,'fro');         
Bww    = G/nrmG;           W   = Bww;                  
EEm    = speye(m); 
nrmFG  = nrmF*nrmG;    

% Initialize the first shift values and compute their LU factorizations
Asnew  = As1;  [LA, UA, PA, QA]  = lu(sparse(A-Asnew*EEn));    
Bsnew  = Bs1;  [LB, UB, PB, QB]  = lu(sparse(B'-Bsnew*EEm));

% --- Main Iteration Loop ---
while iter < 200
    iter                 = iter + 1;
    Avv                  = QA*(UA\(LA\(PA*Avv)));
    Bww                  = QB*(UB\(LB\(PB*Bww)));
    
    % --- Orthonormalization of the Krylov subspace basis vectors ---
    % First pass of classical Gram-Schmidt orthogonalization
    for ii = 1:iter
       Avv_ii            = V(:,(ii-1)*p+1:ii*p);
       HA(ii,iter)       = trace(Avv_ii'*Avv); 
       Avv               = Avv-HA(ii,iter)*Avv_ii;
       
       Bww_ii            = W(:,(ii-1)*p+1:ii*p);
       HB(ii,iter)       = trace(Bww_ii'*Bww); 
       Bww               = Bww-HB(ii,iter)*Bww_ii;
    end
    % Second pass of classical Gram-Schmidt orthogonalization (re-orthogonalization)
    for ii = 1:iter     
        Avv_ii          = V(:,(ii-1)*p+1:ii*p);
        ta              = trace(Avv_ii'*Avv);
        HA(ii,iter)     = HA(ii,iter) + ta; 
        KKAA0(ii,iter)  = HA(ii,iter)*Asnew;
        Avv             = Avv - Avv_ii*ta;
        
        Bww_ii          = W(:,(ii-1)*p+1:ii*p);
        tb              = trace(Bww_ii'*Bww);
        HB(ii,iter)     = HB(ii,iter) + tb; 
        KKBB0(ii,iter)  = HB(ii,iter)*Bsnew;
        Bww             = Bww - Bww_ii*tb;
    end 
    % Update the projected matrix of A(B) on V(W)
    HA(iter+1,iter)     = norm(Avv,'fro');
    Avv                 = Avv/HA(iter+1,iter); 
    V                   = [V,Avv];
    KKAA0(iter+1,iter)  = HA(iter+1,iter)*Asnew;
    AOH                 = [zeros(1,iter-1),HA(iter+1,iter)];
    KKAA1               = KKAA0+[eye(iter);zeros(1,iter)];
    
    HB(iter+1,iter)     = norm(Bww,'fro');
    Bww                 = Bww/HB(iter+1,iter); 
    W                   = [W,Bww];
    KKBB0(iter+1,iter)  = HB(iter+1,iter)*Bsnew;
    BOH                 = [zeros(1,iter-1),HB(iter+1,iter)];
    KKBB1               = KKBB0+[eye(iter);zeros(1,iter)];

    % --- Restart Section (Every km iterations) ---
    if mod(iter,km)==0
        HMA = HA;      HMB = HB; 
        TMA = [];      TMB = [];
        Avm = A*Avv;   Bwm = B'*Bww;
        % First pass of orthogonalization for auxiliary matrix Avm(Bwm)
        for ii = 1:iter+1
           vv_ii           = V(:,(ii-1)*p+1:ii*p);
           HMA(ii,iter+1)  = trace(vv_ii'*Avm); 
           Avm             = Avm-HMA(ii,iter+1)*vv_ii;
           
           ww_ii           = W(:,(ii-1)*p+1:ii*p);
           HMB(ii,iter+1)  = trace(ww_ii'*Bwm); 
           Bwm             = Bwm-HMB(ii,iter+1)*ww_ii;
        end
        % Second pass of orthogonalization for auxiliary matrix (re-orthogonalization) 
        for ii = 1:iter+1    % re-normalize
            vv_ii          = V(:,(ii-1)*p+1:ii*p);
            ta             = trace(vv_ii'*Avm);
            HMA(ii,iter+1) = HMA(ii,iter+1) + ta; 
            TMA            = [TMA; HMA(ii,iter+1)];
            Avm            = Avm - vv_ii*ta;
            
            ww_ii          = W(:,(ii-1)*p+1:ii*p);
            tb             = trace(ww_ii'*Bwm);
            HMB(ii,iter+1) = HMB(ii,iter+1) + tb; 
            TMB            = [TMB; HMB(ii,iter+1)];
            Bwm            = Bwm - ww_ii*tb;
        end 
        hha  = norm(Avm,'fro');
        Avm  = Avm/hha;  
        
        hhb  = norm(Bwm,'fro');
        Bwm  = Bwm/hhb;  
        
    % --- Solve Low-Dimensional Problem ---
    % Construct the projected matrix AKM((BKM) for the low-dimensional problem  
       AKM  = [(KKAA1-TMA*AOH)/HA(1:iter,1:iter),TMA];
       BKM  = [(KKBB1-TMB*BOH)/HB(1:iter,1:iter),TMB];
       fgm  = sparse(iter+1,iter+1); fgm(1,1) = nrmFG;   
    % Solve the low-dimensional Sylvester equation
       Y    = lyap(AKM,BKM',fgm);   

    % Calculate the residual norm using the low-dimensional solution Y
       ttMA = hha*[(-AOH/HA(1:iter,1:iter)),1];
       ttMB = hhb*[(-BOH/HB(1:iter,1:iter)),1];
       res  = sqrt(norm(ttMA*Y,'fro')^2+norm(Y*ttMB','fro')^2);  
    % Calculate the true error norm using the full-form solution Xm 
       Xm   = V(:,1:iter*p+p)*kron(Y,speye(p))*W(:,1:iter*p+p)';  
       err1 = norm(Xm-XX,2);
       E1   = [E1;err1];
    % Calculate the estimated error norm using the low-dimensional solution Y 
       err2 = res/abs(As1+Bs1);
       E2   = [E2;err2];
       
    % Check for convergence
       if res < tol*nrmFG;break;end

    % --- Prepare for Next Cycle (Update Poles and Factorizations) ---
       As2     = norm(AKM,1);    As1  = As2/condest(AKM);
       Bs2     = norm(BKM,1);    Bs1  = Bs2/condest(BKM);
    % Select new poles Asnew and Bsnew based on the current state
       Asnew   = newpole(ttMA,AKM,As1,As2,SSA);
       Bsnew   = newpole(ttMB,BKM,Bs1,Bs2,SSB);
       SSA     = [SSA,Asnew];
       SSB     = [SSB,Bsnew];
    % Update the LU factorizations for the new poles
       [LA, UA, PA, QA]  = lu(sparse(A-Asnew*EEn));
       [LB, UB, PB, QB]  = lu(sparse(B'-Bsnew*EEm));
    end  
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = newpole(tt,K,s1,s2,ss)
% Function to select a new optimal pole within the range [s1, s2]
% that avoids already selected poles in ss

mK = size(K,1);   maxs   = 100; 
e1 = zeros(mK,1); e1(1)  = 1;

% Generate logarithmically spaced test poles between s1 and s2
s_test  = logspace(log10(s1),log10(s2),maxs);

Im      = speye(mK); 
val     = sparse(maxs,1);

% Evaluate the function value at each test pole
for k=1:maxs   
    Ym     = tt*((s_test(k)*Im-K)\e1);  
    val(k) = abs(Ym);
end
[~,id] = sort(val,'descend'); s = abs(s_test(id(1)));

% Check if the selected pole s is already in the list ss (to avoid duplicates)
jk    = 1;
while any (abs(s-ss)<1e-15)  && jk <maxs
jk    = jk+1;
disp('double NGRKSpole found'),
% If duplicate, select the next best pole
s     = abs(s_test(id(jk)));
end
end



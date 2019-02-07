%==========================================================================
% Poisson Solver
% Solves psi_xx + psi_yy = f(x,y) on the x-y grid specified.
% Uses SOR, with the parameter alpha.
% Dirichlet boundary conditions given by psi_init
%==========================================================================
function psi_now = poissonSolverTwoDimensions(f, psi_init, x, y)

%Determine the convergence requirement based on how large the values of psi
%are
convergence_criteria = 10^-7 * (max(psi_init) + 1); %add one incase psi_init is entirely 0
alpha = 1.8;

[dx, dy, x_num, y_num] = meshGridProperties(x,y);

%psi_init=zeros(x_num,y_num);

%Initializing previous and present iterations
psi_now=psi_init;
psi_prev=psi_init;

%Giving initial difference between psi_now and psi_prev to start the
%iterations
psi_prev(round(x_num/2),round(y_num/2))=1;

%Iteration loop
while(max(max(abs(psi_now-psi_prev)))>convergence_criteria)%Run this until convergence
    
    psi_prev=psi_now;

    for x_i=2:(x_num-1)
        
        for y_j=2:(y_num-1)  
            %rhs = gamma * (T_n(x_i, y_j + 1) - T_n(x_i, y_j - 1))/(2*dy);
            psi_now(x_i,y_j)=(1-alpha)*psi_now(x_i,y_j)+(alpha/(2*(dx^2 + dy^2)))*(dy^2 * (psi_now(x_i-1,y_j)+psi_now(x_i+1,y_j))+ dx^2 * (psi_now(x_i,y_j-1)+psi_now(x_i,y_j+1)) - (dx^2*dy^2) * f(x_i,y_j) );
        end
        
    end

end

end
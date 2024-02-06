function calm_convergence_scriptCaller(N, T, lambda, pts)

if ~(exist('N','var')) % Number of x gridpoints
    N = 128;
end
if ~(exist('T','var')) % Final time
    T = 1;
end
if ~(exist('lambda','var')) % Instability parameter
    lambda = 4.1;
end
if ~(exist('pts','var')) % Number of epsilon points
    pts = 11;
end

epsilons = linspace(-14, -4, pts);
epsilons = 10.^(epsilons);

LIL2_error = zeros(6, pts);
L2H2_error = zeros(6, pts);
LILI_error = zeros(6, pts);
   
nametag = ['convergenceData_T_' strrep(num2str(T),'.','p')];
nametag = strrep(nametag,'.','p');

% Test 1: Convergence of type 1 calming with type 1 initial data
for i = 1:pts
   epsilon = epsilons(i);  

   [t1] = KSE2D_IF_RK4_calm(N,T, lambda, epsilon, 1, 1);
   LIL2_error(1,i) = t1.LinfL2_error;
   L2H2_error(1,i) = (t1.L2H2_error_sq)^0.5; % Square root taken because the output is norm squared
   LILI_error(1,i) = t1.LinfLinf_error;
end
disp('Test 1 is done');
save([nametag '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');
disp('Test 1 data is saved');

% Test 2: Convergence of type 2 calming with type 1 initial data
for i = 1:pts
   epsilon = epsilons(i);  

   [t2] = KSE2D_IF_RK4_calm(N,T, lambda, epsilon, 2, 1);
   LIL2_error(2,i) = t2.LinfL2_error;
   L2H2_error(2,i) = (t2.L2H2_error_sq)^0.5; % Square root taken because the output is norm squared
   LILI_error(2,i) = t2.LinfLinf_error;
end
disp('Test 2 is done');
save([nametag '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');
disp('Test 2 data is saved');

% Test 3: Convergence of type 3 calming with type 1 initial data
for i = 1:pts
   epsilon = epsilons(i);  

   [t3] = KSE2D_IF_RK4_calm(N,T, lambda, epsilon, 3, 1);
   LIL2_error(3,i) = t3.LinfL2_error;
   L2H2_error(3,i) = (t3.L2H2_error_sq)^0.5; % Square root taken because the output is norm squared
   LILI_error(3,i) = t3.LinfLinf_error;
end
disp('Test 3 is done');
save([nametag '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');
disp('Test 3 data is saved');

% Test 4: Convergence of type 1 calming with type 2 initial data
for i = 1:pts
   epsilon = epsilons(i);  

   [t4] = KSE2D_IF_RK4_calm(N, T, lambda, epsilon, 1, 2);
   LIL2_error(4,i) = t4.LinfL2_error;
   L2H2_error(4,i) = (t4.L2H2_error_sq)^0.5; % Square root taken because the output is norm squared
   LILI_error(4,i) = t4.LinfLinf_error;
end
disp('Test 4 is done');
save([nametag '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');
disp('Test 4 data is saved');

% Test 5: Convergence of type 2 calming with type 2 initial data
for i = 1:pts
   epsilon = epsilons(i);  

   [t5] = KSE2D_IF_RK4_calm(N,T, lambda, epsilon, 2, 2);
   LIL2_error(5,i) = t5.LinfL2_error;
   L2H2_error(5,i) = (t5.L2H2_error_sq)^0.5; % Square root taken because the output is norm squared
   LILI_error(5,i) = t5.LinfLinf_error;
end
disp('Test 5 is done');
save([nametag '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');
disp('Test 5 data is saved');

% Test 6: Convergence of type 3 calming with type 2 initial data
for i = 1:pts
   epsilon = epsilons(i);  

   [t6] = KSE2D_IF_RK4_calm(N,T, lambda, epsilon, 3, 2);
   LIL2_error(6,i) = t6.LinfL2_error;
   L2H2_error(6,i) = (t6.L2H2_error_sq)^0.5; % Square root taken because the output is norm squared
   LILI_error(6,i) = t6.LinfLinf_error;
end
disp('Test 6 is done');
save([nametag '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');
disp('All test data is saved');


 calm_convergence_plotter(nametag, T, pts);
 
 clear all; 
 close all;
 
end



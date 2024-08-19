%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Medical Imaging Shot Noise Removal 
% ECE565 Term Project
% Poisson noise removal via Expectation Maximization
%
% YongHwan Lee and Tony Storey
%
% given:
%   P(Lambda) = (exp(A_ij*theta_i) * (A_ij*theta_i)^Y_ij)/(Y_ij!)
%   yn = Poisson((A*th)n)   n = 1,2,3...n
%   use EM algorithm    
%   
%
%   
%
%           3x3 Cross section of Voxel model
%           each pixel modeles absorbtion coef
%                    
%    |        *--------------------*       |  
%    | y3---\ |      |      |      |       | 
%    |   ---/ |  p1  |  p2  |  p3  |       | 
%    -        *--------------------*       - 
%    | y2---\ |      |      |      |       | 
%    |   ---/ |  p4  |  p5  |  p6  |       | 
%    -        *--------------------*       - 
%    | y1---\ |      |      |      |       | 
%    |   ---/ |  p7  |  p8  |  p9  |       | 
%    |        *--------------------*       |  
%
%
%             -------|------|-------
%                y9     y10    y11
%                ||     ||     ||
%                \/     \/     \/
%             *--------------------*          
%             |      |      |      |        
%             |  p1  |  p2  |  p3  |        
%             *--------------------*        
%             |      |      |      |        
%             |  p4  |  p5  |  p6  |        
%             *--------------------*        
%             |      |      |      |        
%             |  p7  |  p8  |  p9  |        
%             *--------------------*         
%   
%             -------|------|-------
%

close all
clear
clc

%% Initialize
gain = [0.01 0.1 1 10 100 1000];
max_gain_stages = 6
iter = 30;                                  % Iteration Number
monte = 200;                                 % Montecarlo Number
n = 9;                                      % A matrix length n
m = 16;                                     % A matrix length m

%% Model Matrix A
A = [0 0 0 0 0 0 1 1 1;         
     0 0 0 1 1 1 0 0 0;
     1 1 1 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0;
     0 0 0 1 0 0 0 1 0;
     1 0 0 0 1 0 0 0 1;
     0 1 0 0 0 1 0 0 0;
     0 0 1 0 0 0 0 0 0;
     1 0 0 1 0 0 1 0 0;
     0 1 0 0 1 0 0 1 0;
     0 0 1 0 0 1 0 0 1;
     1 0 0 0 0 0 0 0 0;
     0 1 0 1 0 0 0 0 0;
     0 0 1 0 1 0 1 0 0;
     0 0 0 0 0 1 0 1 0;
     0 0 0 0 0 0 0 0 1];

%% true data + noise for absorption coefs
% Create Observation (Noised)
true_value = [255 270 225 292 218 219 290 235 290].';     % Parameter Vector
x = gain(1)*true_value;
y = A*poissrnd(x);
% x_crlb = inv(A'*A)*A'*y;
% %% Calculating CRLB
% CRLB_init = inv(A.'*diag(1./(A*x_crlb))*A);
% CRLB = diag(CRLB_init);    %./gain(1)
%% Starting point from Method of Moments 
% method of moments starting point
%X = inv(A.'*A)*A.'*y;
    for i = 1:9
      X(i,1) = sum(x)/9;                      % Initial Value
    end
%% EM Algorithm
SE_diag = zeros(9,9);
figure
hold on

p1 = plot (x,'b -o');%true val
p2 = plot (X,'r -s');%init val

% gain value sweep
for gain_val = 1:max_gain_stages
    
    % Monte-Carlo Runs
    for a = 1:monte

        % Observation Update With New Noise
        x = true_value;                         % Parameter Vector
        x = gain(gain_val)*x;
        y = A*poissrnd(x);
        
%         x_crlb = inv(A'*A)*A'*y;
%         CRLB_iter = inv(A.'*diag(1./(A*x_crlb))*A);
 
        
        
        p1 = plot (x,'b -o');%True Value    
        %X = inv(A.'*A)*A.'*y;
        for i = 1:9
          X(i,1) = sum(x)/9;                     % Initial Value
        end
        p2 = plot (X,'r -s'); %init val
        
        % Starting EM Algorithm
        for N = 1:iter                                 
            % M step of EM Algorithm
            for j = 1:n
                for i = 1:m
                    Den1 = A*X;
                    T1(i) = y(i)*A(i,j)/Den1(i);
                end
                Right(j) = sum(T1);
                Den2(j,1) = sum(A(:,j));
                X(j,1) = X(j,1)/Den2(j)*Right(j);
            end

            %Save updated parameters for plotts 
            X_Val(:,N) = X;
            %compute log likelihood 
            likelihood(N) =  sum(y'*log(A*X))-sum(A*X);%-sum(log(factorial(y)));
            
            p3 = plot (X,'y');%iteration      

            % Compute MSE for parameters this poistion sums squared
            % differences
            SE_diag = SE_diag + ((x-X))*((x-X)).';  %diagonal of matrix is summed squared error 
            % Re calculate CRLB with noise 
            CRLB_iter = inv(A.'*diag(1./(A*X))*A);
            % Piece together gain sweep plot info
            SE = diag(SE_diag);                %construct squared error  (:,N)   
            CRLB(:,N) = diag(CRLB_iter);            %construct crlb         
        end% End of EM Algorithm N-iter
        
        
        % find the minimum values for MSE and CRLB per parameter   
        % divide last summed iteration SE by number of iterations to get MSE 
        MSE = SE./iter;  %(:,iter)

        % store an MSE (9x1) vector for each MC iteration
        MSE_monte(:,a) = MSE;
        CRLB_monte(:,a) = CRLB(:,N);   
    
        X_monte(:,a) = X;

    p4 = plot (X,'g -x');%final value                     
    end % end of Monte_Carlo a-iter
   
    % compute the mean normalized values of MSE and CRLB and use
    % to format for log plot normalize by gain factor  
    for par_num = 1:9
        CRLB_monte_min(par_num, gain_val) = log(mean(CRLB_monte(par_num,:)));  %   
        MSE_monte_min(par_num, gain_val) = log(mean(MSE_monte(par_num,:)));  %    
    end                 
 
end % end gain value sweep
%% Plots for convergence, mse, crlb, and log likelihood

hold off
title('Convergence Plot for 9 Voxel Absorption Coefficients')
xlabel('Voxel Absorption Coefficient Number (n)') 
ylabel('Absorption Values Log(p(n))') 
legend([p1, p2, p3, p4], 'True Value', 'Iter = init MM', 'Iter = 1 to 10', 'Final Value', 'Location','southeast')

%plot CRLB after "monte" monte carlo runs for each parameter
figure
for plot_num = 1:9
    subplot(3,3,plot_num)
    hold on
    plot(CRLB_monte_min(plot_num,:), 'b')
    plot(MSE_monte_min(plot_num,:), 'r')    
    hold off
    xlabel(['gain factor ', num2str(gain(1)), '\leq gain \leq', num2str(gain(max_gain_stages))])
    ylabel('MSE & CRLB')
    legend('CRLB', 'MSE')
    title('CRLB&MSE_m_o_n_t_e_m_i_n')
end

figure
plot(likelihood, '-r')
title('Log Likelihood')
xlabel('EM Iteration') 
ylabel('Log Likelihood') 




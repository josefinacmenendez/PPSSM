function firing_rate_struct = get_FR_function(I, J,initial_guess, sigesqguess, sigesqkguess, xkguess)
% INPUTS:
% I: sum of spikes (across channels or neurons)
% J: total number of channels or neurons
% initial_guess: true or false, if true, estimates the initial state by
% flipping the data
% sigesqguess: initial guess for state process variance
% sigesqkguess: initial guess for the one-step prediction variance
% (sig_sq{k|k-1})
% xkguess: initial guess for the one-step prediction mean (x_{k|k-1}) 
% FIXED PARAMETERS (can also be changed)
%       n_iters: number of Newton iterations
%       n_crit: convergence criteria for Newton's method
%       em_steps: maximum number of EM iterations
%       small (which could lead to convergence issues)
%       Parameters specific to the Armijo Line Search:
%           delta: search control parameter (bounded within [0,1])
%           epsilon: search control parameter (bounded within [0,1])
%           alpha: initial Newton step (initialized to 1)
%           alpha_armijo: modified Newton step (found via Backtracking
%           Armijo Line Search)

em_crit = 1e-5;
sigesqguess_init = sigesqguess;
if initial_guess
    I = flip(I);
    S_t = find(I);
    t_s = S_t(end)-1;
    em_crit = 1e-3;
end
%Set parameters for Newton's method---------------
n_iters = 200;
n_crit  = 1e-5; 
delta = 0.5; 
epsilon = 1e-3; 
%--------------------------------------------------------------------------
em_steps  = 200;  %maximum time steps for EM
em_step = 1;
converged=false;
s = 2;             %start recursive filter at trial 2
%sige = sqrt(sig_sq_k_k_guess);% starting guess for state process variance
T                = size(I,2);
x_k_k         = zeros(1,T);
sig_sq_k_k = zeros(1,T);
x_k_k_m1 = zeros(1, T+1);
sig_sq_k_k_m1 = zeros(1, T+1);
%--------------------------------------------------------------------------
% Loop for EM--------------------------------------------------------------
tic;
while em_step < em_steps && ~ converged
%--------------------------------------------------------------------------
% Do forward filter--------------------------------------------------------
x_k_k(s-1)         = xkguess;
sig_sq_k_k(s-1) = sigesqkguess;  
%--------------------------------------------------------------------------
% Loop through all time----------------------------------------------------
for t = s:T+1
   x_k_k_m1(t)      = x_k_k(t-1);
   sig_sq_k_k_m1(t) = sig_sq_k_k(t-1) + sigesqguess;
%--------------------------------------------------------------------------
% Newton's method for finding x_k_k----------------------------------------
    I_t = I(t-1);
    x_k_k_m1_t      = x_k_k_m1(t);
    sig_sq_k_k_m1_t = sig_sq_k_k_m1(t);
    x_k_k_m1_l      = x_k_k_m1_t + (sig_sq_k_k_m1_t*(-exp(x_k_k_m1_t) + I_t));
    f  = @(x_l) 1/(2* sig_sq_k_k_m1_t)*((x_l-x_k_k_m1_t).^2)-I_t*x_l+J*exp(x_l) ;
    g = @(x_l) (x_l-x_k_k_m1_t)./sig_sq_k_k_m1_t-I_t+J*exp(x_l);
    h = @(x_l) 1/sig_sq_k_k_m1_t + J*exp(x_l);
% -------------------------------------------------------------------------
    flagfail=0;
    n_iter  = 1;
    n_cvrgd = false;
    nfail   = 0;
    alpha = 1; alpha_armijo = 1; 
while n_iter < n_iters && ~ n_cvrgd
    d = -g(x_k_k_m1_l)/h(x_k_k_m1_l);
    % Find step size using armijo linesearch method
    j    = 1;
    while (j>0)
        x_new = x_k_k_m1_l+alpha.*d;
        if ( f(x_new) <= f(x_k_k_m1_l)+epsilon*alpha*g(x_k_k_m1_l)*d )
            j = 0;
            alpha_armijo = alpha;
        else
            alpha = alpha*delta;
        end    
    end
   x_k_given_k_l = x_k_k_m1_l + (d*alpha_armijo);
   if abs( (x_k_given_k_l - x_k_k_m1_l)) < n_crit
      %fprintf(2,'Converged (%d) in %d iters \n', n_crit, n_iter);
      flagfail = 0; 
      n_cvrgd = true;
   else
       x_k_k_m1_l = x_k_given_k_l;
       n_iter = n_iter + 1;
   end
end

if n_iter > n_iters && ~ newton_cvrgd
    flagfail = 1;
end
nfail = nfail + flagfail;
x_k_k(t) = x_k_given_k_l;
if nfail > 0
fprintf(2,'Newton convergence failed in %d time points. \n', nfail)
end
%--------------------------------------------------------------------------
% Getting posterior variance-----------------------------------------------
    sig_sq_k_k(t) = 1/( (1/sig_sq_k_k_m1(t)) + (J * exp(x_k_k(t))) );  
end
clear t;
%--------------------------------------------------------------------------
% Do backward filter-------------------------------------------------------
x_k_K        = zeros(1,T);
sig_sq_k_K   = zeros(1,T);
x_k_K(T)     = x_k_k(T);
sig_sq_k_K(T)= sig_sq_k_k(T);
for t = T-1 :-1: s-1
   a(t)         = sig_sq_k_k(t)/sig_sq_k_k_m1(t+1);
   x_k_K(t)     = x_k_k(t) + (a(t)*(x_k_K(t+1) - x_k_k_m1(t+1)));
   sig_sq_k_K(t)= sig_sq_k_k(t)+(a(t)*a(t)*(sig_sq_k_K(t+1)-sig_sq_k_k_m1(t+1)));
end
%--------------------------------------------------------------------------
% Mstep--------------------------------------------------------------------
x_k_K(1) = []; sig_sq_k_K(1) = []; a(1) =[];

N = length(x_k_K);
W_k_K = x_k_K(1:N-1).^2 + sig_sq_k_K(1:N-1);
W_k_plus_one_k_K = ( a.* sig_sq_k_K(2:N) ) + (x_k_K(1:N-1) .* x_k_K(2:N));
W_k_plus_one_K   = x_k_K(2:N).^2 + sig_sq_k_K(2:N);
sigsq_new = sum( W_k_K - (2 *W_k_plus_one_k_K) + W_k_plus_one_K ) / (N-1);
%--------------------------------------------------------------------------
% Set starting values for the next step jk---------------------------------
sigesqguess = sigsq_new;
sigesqkguess = sig_sq_k_K(1);
xkguess = x_k_K(1);
if sigesqguess < 1e-6
    sigesqguess = sigesqguess_init;
end
%--------------------------------------------------------------------------
%check for convergence-----------------------------------------------------
if(em_step>1)
    a1 = abs( (sigsq_new - sigsq_old) )  ;
    if initial_guess
        qnew = x_k_K(t_s); snew = sig_sq_k_K(t_s);
        a2 = abs( (qold - qnew) / qold );
        a3 = abs( (sold - snew)  );
        fprintf('Current absolute changes in the state process variance: %d \n', a1);
        fprintf('Current relative changes in the posterior mode: %d \n', a2);
        fprintf('Current absolute changes in the posterior variance: %d \n', a3);
        converged = (a1 < em_crit) && (a2 < em_crit) && (a3 < em_crit);
        if converged
        firing_rate_struct = struct('spike_counts', I, ...
                            'x_0_K',qnew,...
                            'sig_sq_0_K',snew,...
                            'sige_sq',sigsq_new,...
                            'a', a, 'ccrit', em_crit, ...
                            'converged', converged, ...
                            'em_iters', em_step);
        end
    else
        converged = (a1 < em_crit);
        fprintf('Current absolute changes in the state process variance: %d \n', a1);
        if converged
        firing_rate_struct = struct('spike_counts', I, ...
                            'x_k_given_K',x_k_K,...
                            'sig_sq_k_given_K',sig_sq_k_K,...
                            'sig_epsilon',sigesqguess,...
                            'a', a, 'ccrit', em_crit, ...
                            'converged', converged, ...
                            'em_iters', em_step);
        end
    end
    if( converged )
        t_end = toc;
        fprintf(2,'EM estimates converged in %d iterations and %4.1f seconds \n',em_step, t_end)
    end
end
if initial_guess
    qold = x_k_K(t_s);
    sold = sig_sq_k_K(t_s);
end
sigsq_old = sigsq_new;
%--------------------------------------------------------------------------
em_step = em_step + 1;
end
if em_step == em_steps && a1 > em_crit 
    fprintf(2,'EM failed to converge after %d steps with ccrit %f \n', em_step, em_crit)
    firing_rate_struct = struct('spike_counts', I, ...
                            'x_k_given_K',x_k_K,...
                            'sig_sq_k_given_K',sig_sq_k_K,...
                            'sig_epsilon',sigesqguess,...
                            'a', a, 'ccrit', em_crit, ...
                            'converged', converged, ...
                            'em_iters', em_step);
end

end
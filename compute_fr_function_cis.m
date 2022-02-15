function fr_cis = compute_fr_function_cis(spikes_struct, N_draws, alpha)
% INPUTS: 
%       spikes_struct, the output of 'get_firing_rate'
%       N_draws: number of Monte Carlo samples
%       alpha: confidence level (i.e. 0.05)
x_k_K              = spikes_struct.x_k_given_K;
sig_sq_k_K         = spikes_struct.sig_sq_k_given_K;
a                        = spikes_struct.a;
ucl  = 1 - alpha/2;
lcl   = alpha/2;
lcp  = lcl ;
ucp = ucl;
cps = [lcp 0.5 ucp];
cidx= floor( cps .* N_draws );
L = length(x_k_K);
state_space_covariance_vector = a .* sig_sq_k_K(2:L);
fr_cis = zeros(3,L);
x_lm1 = randn(N_draws,1)* sig_sq_k_K(1) + x_k_K(1);
for l = 1:L-1
    x_lm1 = sort(randn(N_draws,1)*sqrt(sig_sq_k_K(l+1) - ( state_space_covariance_vector(l)^2 / sig_sq_k_K(l+1) )) + x_k_K(l+1) + (( state_space_covariance_vector(l) / sig_sq_k_K(l+1) ) * (x_lm1- x_k_K(l)) )); fr_cis(:,l)=x_lm1(cidx);
end
end
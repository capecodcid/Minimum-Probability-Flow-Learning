% ATH_K_dK_RBM.m
% 
% Author: Andrew Hartnett
% I found K_dK_RBM.m to be difficult to understand
% As I will need to rewrite much of this for different models
% I considered it prudent to make this more understandable
%
% Should write a document to accompany this showing where the
% terms come from.
%
% Additionally this takes a form more akin to model structs we have from our
% other projects
%
% needs to be optimized 2x slower than JS-D version

% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

function [K, dK] = ATH_K_dK_RBM( Wmpf, X )

    [nvis nbatch] = size(X);
    nparams = size(Wmpf,1);
    nhid = nparams/(nvis+1) -1;
    
    Wmpf = reshape(Wmpf, nhid+1, nvis+1);
    W = -Wmpf(1:nhid,1:nvis);
    hbias = -Wmpf(1:nhid,nvis+1)';
    vbias = -Wmpf(nhid+1,1:nvis);
    
    
    
%    nhid = length(model.b);
    
%     W = W';
%     vbias = model.c;
%     hbias = model.b;
    
    K = 0;
    dK = zeros(nhid+1, nvis+1);  % this is so we can fit dK due to bias terms all in the same place
    dK_w_correction_n = zeros(nhid, nvis); 
    dK_w_ij = zeros(nhid,nbatch);  %may need to transpose
    dK_b_j = zeros(nhid,nbatch);
    dK_a_n = zeros(1,nvis);
    
    term_wb_x        = 1 + exp(W*X + repmat(hbias',1,nbatch));
    d_log_term_wb_x  = exp(W*X + repmat(hbias',1,nbatch)) ./ term_wb_x;
    
for d = 1:nvis
    term_wb_xprime  = 1 + exp(W*X + W(:,d)*(1-2*X(d,:)) + repmat(hbias',1,nbatch));
    d_log_term_wb_xprime = exp(W*X + W(:,d)*(1-2*X(d,:)) + repmat(hbias',1,nbatch))./ term_wb_xprime;
    
    energy_diff_wb = log( prod( term_wb_xprime ./ term_wb_x , 1 ) );
    energy_diff_a = vbias(d)*(1-2*X(d,:));
    K_per_data_point = exp( 0.5 * (energy_diff_a + energy_diff_wb));
    
    K = K + sum( K_per_data_point ); 
    
    dK_w_ij = dK_w_ij + (ones(nhid,1) * K_per_data_point) .* (d_log_term_wb_xprime - d_log_term_wb_x);  % good - exact match
    dK_w_correction_n(:,d) = d_log_term_wb_xprime * (K_per_data_point .* (1 - 2*X(d,:)))';  % good - exact match 
    dK_b_j =  dK_b_j + (ones(nhid,1) * K_per_data_point) .* (d_log_term_wb_xprime - d_log_term_wb_x);  % don't need this 
    dK_a_n(d) = 0.5 * K_per_data_point * (1-2*X(d,:))';   % good - exact match
end

dK(1:nhid,1:nvis) = dK_w_ij * X' + dK_w_correction_n;  % putting the dk/dw terms in the matrix
dK(1:nhid,nvis+1) = dK_b_j * ones(nbatch,1);
dK = dK/2;
dK(nhid+1,1:nvis) = dK_a_n;


K = K/nbatch;
dK = -dK / nbatch;

dK = dK(:);

end
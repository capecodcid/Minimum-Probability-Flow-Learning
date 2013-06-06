% - trains an RBM using MPF (code modified from JS-D)
% - original licence below
%
% Andrew Hartnett (2013) ahartnet@princeton.edu
%

% Author: Jascha Sohl-Dickstein (2010)
% Web: http://redwood.berkeley.edu/wiki/Jascha_Sohl-Dickstein
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

function J_gamma = ising_gamma_MPF( X, n_gamma, n_learns)

addpath ../3rd_party_code/minFunc

[ndim, nbatch] = size( X );

fprintf( 'initializing training weights\n' );

J = randn( ndim, ndim ) / sqrt(ndim) * 3.;
J = J + J';                             % !ath - symmetrize J
J = J/2;
J = J - diag(diag(J)); % set the diagonal so all the units are 0 bias
J = J - diag(sum(J));
gamma = rand(n_gamma,1) / sqrt(ndim);
J_gamma = [J(:);gamma];

trip = nchoosek((1:ndim),3);                   
gamma_all = trip(threepointkind(trip,5,2)==2,:); 


fprintf( 'estimating weight matrix via MPF\n' );
% MPF
%Wmpf = Winit;
minf_options = [];
maxlinesearch = 10;
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;
%Wmpf = Wmpf(:);

if n_learns < 1
    error('error - need integer number of learns')
end

%minf_options.maxiter = 1000;
for i = 1:n_learns
    learn_t = tic();
    J_gamma = minFunc( @ATH_K_dK_ising_gamma, J_gamma, minf_options, X , gamma_all);
    learn_t = toc(learn_t)
    
    fprintf('learning round %i of %i took %f seconds \n',i ,n_learns, learn_t)
end

%Wmpf = reshape(Wmpf, d_hid+1, d_vis+1);
%model.W = -Wmpf(1:d_hid,1:d_vis)';
%model.b = -Wmpf(1:d_hid,d_vis+1)';
%model.c = -Wmpf(d_hid+1,1:d_vis);
end
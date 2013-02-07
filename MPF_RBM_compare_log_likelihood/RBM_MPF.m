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

function [model W_mpf] = RBM_MPF(d_vis, d_hid, X, n_learns)

addpath ../3rd_party_code/minFunc

%d_vis = 24; % number of units in the visible layer
%d_hid = 20; % number of units in the hidden layer


% batch_size = 100; % number of training samples to generate
% 
% fprintf( 'choosing random RBM\n' );
% tic();
% Wtrue = 4*randn( d_hid+1, d_vis+1 ) / sqrt(d_vis+1);
% independent_steps = 200; % how many steps to go between samples
% fprintf( 'generating data samples\n' );
% X = sample_RBM( Wtrue, batch_size, independent_steps, independent_steps, rand( d_vis, 1 ) > 0.5 );
% X = X > rand(size(X));

%X = bint;  % I use my raw data rather than their randomly generated data

%n_learns = 1;

fprintf( 'initializing training weights\n' );
Winit = randn( d_hid+1, d_vis+1 ) / sqrt(d_vis+1);
%Winit = Wmpf;


fprintf( 'estimating weight matrix via MPF\n' );
% MPF
Wmpf = Winit;
minf_options = [];
maxlinesearch = 100;
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;
Wmpf = Wmpf(:);

if n_learns < 1
    error('error - need integer number of learns')
end

%minf_options.maxiter = 1000;
for i = 1:n_learns
    learn_t = tic();
    Wmpfn = minFunc( @K_dK_RBM, Wmpf(:), minf_options, X );
    Wmpf(:) = Wmpfn(:);
    W_mpf{i} = Wmpf(:);
    learn_t = toc(learn_t)
    
    fprintf('learning round %i of %i took %f seconds \n',i ,n_learns, learn_t)
end

Wmpf = reshape(Wmpf, d_hid+1, d_vis+1);
model.W = -Wmpf(1:d_hid,1:d_vis)';
model.b = -Wmpf(1:d_hid,d_vis+1)';
model.c = -Wmpf(d_hid+1,1:d_vis);
end
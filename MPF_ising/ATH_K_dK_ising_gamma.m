% ATH_K_dK_ising_gamma.m
% Author: Andrew Hartnett
%
% This computes the objective function for an ising model + a small
% number of gamma (third order) terms that are important for our fish
% data.
%
%
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)


% for now let us act as if we only have one 48 bit data vector
% so we are going to seperate the deltaE term into two terms -
% the usual ising term and the gamma term



function [K dK] = ATH_K_dK_ising_gamma( J_gamma, X , gamma_all )

% general
[ndim, nbatch] = size( X );
%trip = nchoosek((1:ndim),3);
%gamma_all = trip(threepointkind(trip,5,2)==2,:);
gamma = zeros(length(gamma_all),1);

J_gamma = J_gamma(:);
J = J_gamma(1:ndim*ndim);
gamma = J_gamma((ndim*ndim)+1:end);

%J = J_gamma;
J = reshape( J, [ndim, ndim] );
J = (J + J')/2;
h = diag(J);
J = J - diag(diag(J));

Ktot = 0;
dKtot = zeros(2448,1);    % replace with length of J_gamma

for k = 1:nbatch;
    x=X(:,k);
    
    % for finding the relavent gammas
    bits_on = find(x');
    E_gamma_x = 0;
    if length(bits_on) >= 3
        data_trip = nchoosek(bits_on,3);
        [~, ~, gamma_inds] = intersect(data_trip,gamma_all,'rows');
        E_gamma_x = -sum(gamma(gamma_inds));
    end
    

    % now we need to find gammas for each neighboring state
    deltaE_gamma = zeros(ndim,1);
    deltaE_ising = zeros(ndim,1);

    %ising part
    for n = 1:ndim
        j_temp = 0;
        for i = 1:ndim
            j_temp = j_temp + J(i,n)*x(i);
        end
        j_temp = j_temp - J(n,n)*x(n);
        deltaE_ising(n) = (2*j_temp + h(n))*(1-2*x(n));
    end

    %gamma part
    if length(bits_on) >= 2
        for n = 1:ndim
            x_ = x;
            x_(n) = abs(x(n)-1);
            bits_on = find(x_');
            E_gamma_x_ = 0;
            if length(bits_on) >=3
                data_trip = nchoosek(bits_on,3);
                [~, ~, gamma_inds] = intersect(data_trip,gamma_all,'rows');
                E_gamma_x_ = -sum(gamma(gamma_inds));  % -0.0408
            end
            deltaE_gamma(n) = E_gamma_x - E_gamma_x_;
        end
    end
    deltaE = deltaE_ising + deltaE_gamma;
    Kn = exp(-0.5*deltaE);
    K = sum(Kn);

    % for the derivative - only bit flips corresponding to a given gamma will
    % count
    % i can loop through the gammas?

    d_deltaE_gamma = zeros(ndim,length(gamma_all));
    dK_J = zeros(ndim, ndim);

    for n = 1:ndim
        temp = zeros(ndim,1);
        j_temp = 0;
        for i = 1:ndim;
            j_temp = j_temp + J(i,n)*x(i);
            temp(i) = (-0.5 * (1-2*x(n)) * 2 * x(i));
        end
        j_temp = j_temp - J(n,n)*x(n);
        termA = (2*j_temp + h(n))*(1-2*x(n));
        temp(n) = (-0.5 * (1-2*x(n)));
        dK_J(:,n) = dK_J(:,n) + exp(-0.5*termA)*temp;
    end

    for g = 1:length(gamma_all)
        if sum(x(gamma_all(g,:)))==2
            off_bit = gamma_all(g,~x(gamma_all(g,:)));
            d_deltaE_gamma(off_bit,g) = -(0.5)*gamma(g);
            
        elseif sum(x(gamma_all(g,:)))==3
            on_x = 1;
            for n = 1:ndim
                on_x_ = 1;
                if ismember(n,gamma_all(g,:)), on_x_ = 0; end
                d_deltaE_gamma(n,g) = -(0.5)*gamma(g)*(on_x_ - on_x);
            end
        end
    end
    dK_J = (dK_J + dK_J') / 2;
    dK_gamma = repmat(Kn,1,length(gamma_all)).*d_deltaE_gamma;
    dK_gamma = sum(dK_gamma,1)';
    dK_J = dK_J(:);
    dK = [dK_J; dK_gamma];
    Ktot = Ktot+K;
    dKtot = dKtot +dK;
    
end
K = Ktot/nbatch;
dK = dKtot/nbatch;

end

% may want to prefilter data - to only run through loops with those with 2+
% activated ( or better yet at least one in each half)


% at bare minimum we can check the ising part of the code by submitting a
% a model with gammas = 0 and comparing K and dK to those computed by the
% regular ising function
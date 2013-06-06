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



function [K dK] = ATH_K_dK_ising_gamma_vector( J_gamma, X )

% general
[ndim, nbatch] = size( X );
trip = nchoosek((1:ndim),3);
gamma_all = trip(threepointkind(trip,5,2)==2,:);            % perhaps just hard code these guys
gamma = 0.01*zeros(length(gamma_all),1);

J_gamma = J_gamma(:);
J = J_gamma(1:(ndim*ndim));
gamma = J_gamma(ndim*ndim+1:end);
gamma = repmat(gamma,1,ndim);
size(gamma)

%J = J_gamma;
J = reshape( J, [ndim, ndim] );
J = (J + J')/2;
h = diag(J);
J = J - diag(diag(J));

dKtot = zeros(2448,1);


% for finding the relavent gammas
bits_on = cellfun(@find,num2cell(X,1),'UniformOutput',0);
num_on = cellfun(@length,bits_on);

E_gamma_x_k(num_on >=3) = cellfun(@(x) -sum(gamma(sum(ismember(gamma_all,x),2)==3),1), bits_on(num_on >=3));

deltaE_ising = (1-2*X).*(2*J*X + repmat(h,1,nbatch));

E_gamma_x_prime = zeros(ndim,nbatch);
t1 = tic();

for n = 1:ndim
    canidate_x_prime_turnoff = cellfun(@(x) ismember(n,x), bits_on);
    %canidate_x_prime_turnoff = cell2mat(cellfun(@(x) ismember((1:ndim),x), bits_on));
    
    x_prime = cellfun(@(x) x(x~=n), bits_on, 'UniformOutput',0);
    x_prime(~canidate_x_prime_turnoff) = cellfun(@(x) sort([x;n]), bits_on(~canidate_x_prime_turnoff), 'UniformOutput', 0);
    %x_prime_on(~canidate_x_prime_turnoff) = cellfun(@(x) sort([x;n]), bits_on(~canidate_x_prime_turnoff), 'UniformOutput', 0);
    %x_prime_off(canidate_x_prime_turnoff) = cellfun(@(x) x(x~=n), bits_on(canidate_x_prime_turnoff), 'UniformOutput', 0);
    
    %num_on_x_prime_on = cellfun(@length,x_prime_on);
    %x_prime_off(canidate_x_prime_turnoff) = cellfun(@(x) x(x~=n), bits_on(canidate_x_prime_turnoff), 'UniformOutput', 0);
    %num_on_x_prime_off = cellfun(@length,x_prime_off);
    num_on_x_prime = cellfun(@length,x_prime);
    
   
    
    %E_x_prime_on = cellfun(@(x) -sum(gamma(sum(ismember(gamma_all,x),2)==3)), x_prime_on(num_on_x_prime_on >=3));
    %E_x_prime_off = cellfun(@(x) -sum(gamma(sum(ismember(gamma_all,x),2)==3)), x_prime_off(num_on_x_prime_off >=3));
    E_gamma_x_prime(n,num_on_x_prime >=3) = cellfun(@(x) -sum(gamma(sum(ismember(gamma_all,x),2)==3,n)), x_prime(num_on_x_prime >=3));
    
    
end
t1 = toc(t1)


E_gamma_x_prime = cell(ndim,1)

t4 = tic();

for n = 1:ndim
    %E_gamma_x_prime{n
    canidate_x_prime_turnoff = cellfun(@(x) ismember(n,x), bits_on);
    
    
    x_prime = cellfun(@(x) x(x~=n), bits_on, 'UniformOutput',0);
    x_prime(~canidate_x_prime_turnoff) = cellfun(@(x) sort([x;n]), bits_on(~canidate_x_prime_turnoff), 'UniformOutput', 0);
    
    num_on_x_prime = cellfun(@length,x_prime);
    
    E_gamma_x_prime(n,num_on_x_prime >=3) = cellfun(@(x) -sum(gamma(sum(ismember(gamma_all,x),2)==3,n)), x_prime(num_on_x_prime >=3));
    
end
t4 = toc(t4)






t3 = tic();
%E = zeros(ndim,nbatch)
E = cell(nbatch,1);
parfor k = 1:nbatch
   E{k} = zeros(ndim,1);
   x = X(:,k);
   for n = 1:ndim
       x_prime = [];
       bits_on = find(x);
       num_on = length(bits_on);
       if ismember(n, bits_on);
        x_prime = bits_on(bits_on~=n);
        num_on = num_on-1;
       else
        x_prime = sort([x_prime;n]);
        num_on = num_on + 1;
       end
       E{k}(n) = -sum(gamma(sum(ismember(gamma_all,x_prime))==3));
       
    
   end
end
t3 = toc(t3)

deltaE_gamma = repmat(E_gamma_x_k,ndim,1) - E_gamma_x_prime;

deltaE = deltaE_ising + deltaE_gamma;
Knx = exp(-0.5*deltaE);
Kn = sum(Knx,2);
K = sum(Kn);


t2 = tic();
parfor k = 1:nbatch
% for the derivative - only bit flips corresponding to a given gamma will
% count
% i can loop through the gammas?
x = X(:,k);

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
        d_deltaE_gamma(off_bit,g) = -(0.5)*gamma(g,1);
        
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
%Ktot = Ktot+K;
dKtot = dKtot +dK;

end

t2 = toc(t2)
K = K/nbatch;
dK = dKtot/nbatch;

end

% may want to prefilter data - to only run through loops with those with 2+
% activated ( or better yet at least one in each half)


% at bare minimum we can check the ising part of the code by submitting a
% a model with gammas = 0 and comparing K and dK to those computed by the
% regular ising function
X = Xall;
% ATH_K_dK_ising.m
% Author: Andrew Hartnett
%
% this is an alternative to the function by J.S-D. and U.K.
% it was a learning exercise for me - and I show several different
% ways of calculating the same quantities.
%
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

%
% General stuff
%

[ndim, nbatch] = size( X );

Jt = reshape( J_0, [ndim, ndim] );
Jt = (Jt + Jt')/2;

h = diag(Jt);
Jt = Jt - diag(diag(Jt));

%
% Vectorized method #1
%

deltaE = (1-2*X).*(2*Jt*X + repmat(h,1,nbatch));
Kns = exp(-0.5*deltaE);
K = sum(Kns(:));

dJ = 2*X';
dh = 1;
dK = Kns.*(-0.5*(1-2*X));
dKJ = dK*dJ;
dKh = diag(sum(-0.5*Kns,2));

dKtotal = dKJ + dKh;
dKtotal = (dKtotal + dKtotal')/2;

dK = dKtotal/nbatch
K = K/nbatch


%
% Vectorized method #2
%

deltaE = (1-2*X).*(2*Jt*X + repmat(h,1,nbatch));
Kns = exp(-0.5*deltaE);
K = sum(Kns(:));                % same as method 1


dJ2 = Kns.*(-0.5*(1-2*X))*2*X';
dH2 = diag(sum(Kns.*(-0.5*(1-2*X)),2));

dt1 = dKJ-diag(diag(dKJ)) + dH2;
dt1 = (dt1+dt1')/2;
dK = dt1/nbatch;
K = K/nbatch;

% %
% % Semi-vectorized method #3  - not finished/functioning
% %
% 
% K = 0;
% dK = zeros(ndim,ndim);
% dK_dE_h = zeros(ndim,1);
% dK_dE_J = zeros(ndim,ndim);
% dK_h = zeros(ndim,1);
% dK_J = zeros(ndim,ndim);
% 
% 
% for d = 1:ndim
%     deltaE = (1-2*X(d,:)).*(repmat(h(d),1,nbatch) + 2*Jt(d,:)*X);
%     Kn = sum(exp(-0.5*deltaE));
%     K = K + Kn;
%     
%     dK_dE_h(d) = sum((1-2*X(d,:)));
%     dK_dE_J(:,d) = sum(2*ones(ndim,1)*(1-2*X(d,:)).*X,2);
%     
%     dK_h(d) = 0.5*Kn* dK_dE_h(d);
%     dK_J(:,d) = 0.5 * Kn * dK_dE_J(:,d);
% end


%
% Explicit loops method #4
%


Obj = 0;                    % K loop

for k = 1:nbatch
    x = X(:,k);
    for d = 1:ndim
        j_temp = 0;
        for i = 1:ndim
            j_temp = j_temp + Jt(i,d)*x(i);
        end
        j_temp = j_temp - Jt(d,d)*x(d);
        termA = (2*j_temp + h(d))*(1-2*x(d));
        
        Obj = Obj + exp(-0.5*termA);
    end
end

dloop = zeros(ndim,ndim);        % dK loop

for k = 1:nbatch
    x = X(:,k);
    for d = 1:ndim
        temp = zeros(ndim,1);
        j_temp = 0;
        for i = 1:ndim;
            j_temp = j_temp + Jt(i,d)*x(i);
            temp(i) = (-0.5 * (1-2*x(d)) * 2 * x(i));
        end
        j_temp = j_temp - Jt(d,d)*x(d);
        termA = (2*j_temp + h(d))*(1-2*x(d));
        temp(d) = (-0.5 * (1-2*x(d)));
        dloop(:,d) = dloop(:,d) + exp(-0.5*termA)*temp;
    end
end

dK = (dloop + dloop')/(2*nbatch);
K = Obj/nbatch;         

%
% Super written out - want to do loop over i and loop over J so that 
%

dloop = zeros(ndim,ndim);        % dK loop

for k = 1:nbatch
    x = X(:,k);
    for d = 1:ndim
        temp = zeros(ndim,1);
        j_temp = 0;
        for i = 1:ndim;
            j_temp = j_temp + Jt(i,d)*x(i);
            temp(i) = (-0.5 * (1-2*x(d)) * 2 * x(i));
        end
        j_temp = j_temp - Jt(d,d)*x(d);
        termA = (2*j_temp + h(d))*(1-2*x(d));
        temp(d) = (-0.5 * (1-2*x(d)));
        dloop(:,d) = dloop(:,d) + exp(-0.5*termA)*temp;
    end
end




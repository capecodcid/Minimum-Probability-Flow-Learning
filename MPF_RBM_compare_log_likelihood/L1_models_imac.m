
nhid = 20;
num_learns = 50;
lambdas = [0.00 1e-5 1e-4 1e-3];

d_vis =24;
data = shuffle_bint(bint);
fold = 4;

data_folds = cell(fold,1);
for i = 1:fold
    data_folds{i} = data(:,i:fold:end);
end

for i = 1:fold
    train = [];
    for ict = 1:fold
        if ict == i
            xval = data_folds{ict};
        else
            train = [train, data_folds{ict}];
        end
    end

    for hct = 1:length(lambdas)
        lambda = lambdas(hct);
        Winit = randn( nhid+1, d_vis+1 ) / sqrt(d_vis+1);
        
        for nl = 1:num_learns
            progress_string = ['Learning round ',num2str((i-1)*(length(lambdas)*num_learns)+(hct-1)*num_learns + nl),' of ',num2str(fold*length(lambdas)*num_learns),' total rounds'];
            disp(progress_string)
            [model Winit] = RBM_MPF_mbp_L1(d_vis, nhid, Winit, train, lambda);
        end
        
        savename = ['lambda_',num2str(hct),'_nhid_',num2str(nhid),'_fold_', num2str(i),'_nl_',num2str(num_learns),'.mat'];
        save(savename,'model','Winit','train','xval')
       

    end


end
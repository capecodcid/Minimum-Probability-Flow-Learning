% sandbox.m


%% trying to figure out which method is faster for finding the id's in each row

t1 = tic();

for i = 1:100
C = num2cell(bint,1);             %# Collect the columns into cells
ids = cellfun(@find,C,'UniformOutput',0);
end
t1 = toc(t1)

applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :))
applyToRows = @(func, matrix) arrayfun(applyToGivenRow(func, matrix), 1:size(matrix,1), 'UniformOutput',0)'


t2 = tic();
for i = 1:100
x = applyToRows(@find,bint');
end
t2 = toc(t2)

%%

for k = 1:1000
    X_ = bint(:,k);
    [K_ , dK_] = ATH_K_dK_ising_gamma( J_0, X_);
    Kath = Kath + K_;
    dKath = dKath + dK_;
end
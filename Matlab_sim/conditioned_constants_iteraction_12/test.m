% para
ll = zeros(10,1);

for i = 1:10
    para = [i; para_est(2:10)];
    ll(i) = log_ll(para, wait_cost, beta, T, N, ...
                    alloc_vec, category, tran_matrix, time_index, state_index, cross_clamp);

end
plot([1:10], ll)
function obj=log_ll(para, wait_cost, beta, T, N, ...
                    alloc_vec, category, tran_matrix, time_index, state_index, cross_clamp)
%% Compute choice prob for all states and all t                

% para = para_est(1:6);
para_vec=cell(4,1);
para_vec{1,1}=para(1:6);
para_vec{2,1}=[para(1); para(3); para(5)];
para_vec{3,1}=[para(2); para(3); para(6)];
para_vec{4,1}=para(3);

V_vec=cell(4,1);
choice_prob_vec=cell(4,1);

temp = alloc_vec{1,1};
interaction = [temp(:,1).*temp(:,2) temp(:,1).*temp(:,3) temp(:,2).*temp(:,3)];
alloc_vec{1,1} = [temp interaction];

temp = alloc_vec{2,1};
interaction = [temp(:,1).*temp(:,2)];
alloc_vec{2,1} = [temp interaction];

temp = alloc_vec{3,1};
interaction = [temp(:,1).*temp(:,2)];
alloc_vec{3,1} = [temp interaction];
costs = [-1; para(7)];
c = para(8);

for j=1:2
%     c = para(6+j);
    wait_cost = costs(j);
    para_temp=para_vec{j,1};
    alloc_vec_temp=alloc_vec{j,1};
%     tran_matrix_temp=tran_matrix{j,1};    
    n_comb=size(alloc_vec_temp,1);
    
    
    V_temp=zeros(T,n_comb);
    choice_prob_temp=zeros(T,n_comb);
    
    % flow utility if chooses cross clamp
    v1=(alloc_vec_temp*para_temp)';
    
    
    v1 = v1+c;
    
    % at T, always chooses cross clamp
    V_temp(T,:)=v1;
    

    for t=T-1:-1:1
        tran_matrix_temp=tran_matrix{j,t};
        v0=wait_cost+beta*sum(repmat(V_temp(t+1,:), [n_comb,1]).*tran_matrix_temp,2)';    
        V_temp(t, :)=0.5772+log(exp(v0)+exp(v1));
        choice_prob_temp(t,:)=exp(v0)./(exp(v0)+exp(v1));
        
    end
    
    V_vec{j,1}=V_temp;
    choice_prob_vec{j,1}=choice_prob_temp;
    
end
                
%% Construct likelihood for choices

p_vec=zeros(N,2);
for i=1:N
    
    choice_prob_temp=choice_prob_vec{category(i),1};
    p0=choice_prob_temp(time_index(i),state_index(i));
    p_vec(i,:)=[p0, 1-p0];
    
end

obj=-sum(log(prod(p_vec.^[1-cross_clamp,cross_clamp],2)),1)/N;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal Stopping Problem
% 05/03/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%


% clc
% clear
% 
% % path='/Users/johnwang/OneLegacy/code_Yi/matlab_code';
% % cd(path);
% 
% 
% load('est_data_1h_v2.mat');



% time interval between two obs: 1 hour
% ignore intestine to reduce state space
% when cleaning the dataset: one donor has LU alone, dropped


%% load data and set parameters

organ_available=[HR_available, LU_available, LI_available];   
organ_allocated=[HR_allocated, LU_allocated, LI_allocated];


wait_cost=-1;
beta=0.99;
T=max(time_index);  % last period
n=max(donor_index); % number of donors
n_org=size(organ_available, 2); % number of organs considered

N=size(donor_index, 1); % number of rows


alloc_vec_c1=combvec([0,1], [0, 1, 2], [0,1])';  % category=1: HR+LU+LI
alloc_vec_c2=combvec([0,1], [0,1])';             % category=2: HR+LI
alloc_vec_c3=combvec([0, 1, 2], [0,1])';         % category=3: LU+LI
alloc_vec_c4=combvec([0,1])';                    % category=4: LI

HR_vec=combvec([0,1], [0,1])';
LU_vec=combvec([0,1,2], [0,1,2])';
LI_vec=combvec([0,1], [0,1])';

%% first estimate Pr(organ is allocated at t+1 | organ is not allocated at t)
% do this for each organ

% % HR
% temp=zeros(N,2);
% for i=1:N-1    
%     temp(i, 1)=(organ_available(i, 1)==1 & organ_allocated(i, 1)==0 & donor_index(i)==donor_index(i+1));
%     temp(i, 2)=(organ_available(i, 1)==1 & organ_allocated(i, 1)==0 & organ_allocated(i+1, 1)==1 & donor_index(i)==donor_index(i+1));
% end
% temp=sum(temp, 1);
% prob_HR=[1-temp(2)/temp(1); 0; temp(2)/temp(1); 1];
% 
% 
% % LU
% temp=zeros(N,5);
% for i=1:N-1    
%     temp(i, 1)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==0 & donor_index(i)==donor_index(i+1));
%     temp(i, 2)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==0 & organ_allocated(i+1, 2)==1 & donor_index(i)==donor_index(i+1));
%     temp(i, 3)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==0 & organ_allocated(i+1, 2)==2 & donor_index(i)==donor_index(i+1));
%     
%     temp(i, 4)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==1 & donor_index(i)==donor_index(i+1));
%     temp(i, 5)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==1 & organ_allocated(i+1, 2)==2 & donor_index(i)==donor_index(i+1));
%     
% end
% temp=sum(temp, 1);
% prob_LU=[1-(temp(2)+temp(3))/temp(1); 0; 0; temp(2)/temp(1); 1-temp(5)/temp(4); 0; temp(3)/temp(1); temp(5)/temp(4); 1];
% 
% 
% 
% % LI
% temp=zeros(N,2);
% for i=1:N-1    
%     temp(i, 1)=(organ_available(i, 3)==1 & organ_allocated(i, 3)==0 & donor_index(i)==donor_index(i+1));
%     temp(i, 2)=(organ_available(i, 3)==1 & organ_allocated(i, 3)==0 & organ_allocated(i+1, 3)==1 & donor_index(i)==donor_index(i+1));
% end
% temp=sum(temp, 1);
% prob_LI=[1-temp(2)/temp(1); 0; temp(2)/temp(1); 1];




%% Construct transition matrix

tran_matrix=cell(4,T);
HR_param = [-2.19 -0.0546];
% LU_param = [-3.51 -0.0085];
% LI_param = [-3.04 -0.0094];
    
for k = 1:T
    
    p_HR = exp(HR_param(1) + HR_param(2) * k)/(1+exp(HR_param(1) + HR_param(2) * k));
%     p_LU = exp(LU_param(1) + LU_param(2) * k)/(1+exp(LU_param(1) + LU_param(2) * k));
%     p_LI = exp(LI_param(1) + LI_param(2) * k)/(1+exp(LI_param(1) + LI_param(2) * k));
    p_LU = 0.042;
    p_LI = 0.0181;
    prob_HR = [1-p_HR p_HR];
    prob_LU = [(1-p_LU) 0 p_LU];
    prob_LI = [(1-p_LI) p_LI];
    prob_LU_0 = [(1-0.012-0.0055) 0.0055 0.012 (1-0.0096) 0.0096];
    prob_LU_1 = [(1-0.037-0.0046) 0.0046 0.037 (1-0.0204) 0.0204];
    prob_LI_0 = [(1-0.0182) 0.0182];
    prob_LI_1 = [(1-0.0759) 0.0759];
    
    tran_matrix_c1=zeros(size(alloc_vec_c1,1), size(alloc_vec_c1,1));

    for i=1:size(alloc_vec_c1,1)
        for j=1:size(alloc_vec_c1,1)
            a = alloc_vec_c1(j,:);
            b = alloc_vec_c1(i,:);
            diff = a - b;
            if ~all(diff >= 0)
                tran_matrix_c1(i, j)=0;
            else
                if b(1) == 0
                    if b(2) == 0
                        if b(3) == 0
                            tran_matrix_c1(i, j)=prob_HR(diff(1)+1)*prob_LU_0(diff(2)+1)*prob_LI_0(diff(3)+1);
                        else 
                            tran_matrix_c1(i, j)=prob_HR(diff(1)+1)*prob_LU_0(diff(2)+1);
                        end
                    elseif b(2) == 1
                        if b(3) == 0
                            tran_matrix_c1(i, j)=prob_HR(diff(1)+1)*prob_LU_0(diff(2)+4)*prob_LI_0(diff(3)+1);
                        else 
                            tran_matrix_c1(i, j)=prob_HR(diff(1)+1)*prob_LU_0(diff(2)+4);
                        end
                    else
                        if b(3) == 0
                            tran_matrix_c1(i, j)=prob_HR(diff(1)+1)*prob_LI_0(diff(3)+1);
                        else 
                            tran_matrix_c1(i, j)=prob_HR(diff(1)+1);
                        end
                    end
                else 
                    if b(2) == 0
                        if b(3) == 0
                            tran_matrix_c1(i, j)=prob_LU_1(diff(2)+1)*prob_LI_1(diff(3)+1);
                        else 
                            tran_matrix_c1(i, j)=prob_LU_1(diff(2)+1);
                        end
                    elseif b(2) == 1
                        if b(3) == 0
                            tran_matrix_c1(i, j)=prob_LU_1(diff(2)+4)*prob_LI_1(diff(3)+1);
                        else 
                            tran_matrix_c1(i, j)=prob_LU_1(diff(2)+4);
                        end
                    else
                        if b(3) == 0
                            tran_matrix_c1(i, j)=prob_LI_1(diff(3)+1);
                        else 
                            tran_matrix_c1(i, j)=1;
                        end
                    end
                end
                
            end
            
        end
    end


    tran_matrix_c2=zeros(size(alloc_vec_c2,1), size(alloc_vec_c2,1));

    for i=1:size(alloc_vec_c2,1)
        for j=1:size(alloc_vec_c2,1)
            a = alloc_vec_c2(j,:);
            b = alloc_vec_c2(i,:);
            diff = a - b;
            if ~all(diff >= 0)
                tran_matrix_c2(i, j)=0;
            else
                if b(1) == 0
                    if b(2) == 0
                        tran_matrix_c2(i, j)=prob_HR(diff(1)+1)*prob_LI_0(diff(2)+1);
                    else
                        tran_matrix_c2(i, j)=prob_HR(diff(1)+1);
                    end
                else
                    if b(2) == 0
                        tran_matrix_c2(i, j)=prob_LI_1(diff(2)+1);
                    else
                        tran_matrix_c2(i, j)=1;
                    end
                end
            end
        end
    end


    tran_matrix_c3=zeros(size(alloc_vec_c3,1), size(alloc_vec_c3,1));

    for i=1:size(alloc_vec_c3,1)
        for j=1:size(alloc_vec_c3,1)
            a = alloc_vec_c3(j,:);
            b = alloc_vec_c3(i,:);
            diff = a - b;
            if ~all(diff >= 0)
                tran_matrix_c3(i, j)=0;
            else
                if b(1) == 0
                    if b(2) == 0
                        tran_matrix_c3(i, j)=prob_LU(diff(1)+1)*prob_LI(diff(2)+1);
                    else 
                        tran_matrix_c3(i, j)=prob_LU(diff(1)+1);
                    end
                elseif b(1) == 1
                    if all(diff == 0)
                        tran_matrix_c3(i,j)=1;
                    else
                        tran_matrix_c3(i, j)=0;
                    end
%                     else 
%                         tran_matrix_c1(i, j)=prob_LU(diff(1)+4);
%                     end
                else
                    if b(2) == 0
                        tran_matrix_c3(i, j)=prob_LI(diff(2)+1);
                    else 
                        tran_matrix_c3(i, j)=1;
                    end
                end
            end
        end
    end


    tran_matrix_c4=zeros(size(alloc_vec_c4,1), size(alloc_vec_c4,1));

    for i=1:size(alloc_vec_c4,1)
        for j=1:size(alloc_vec_c4,1)
            a = alloc_vec_c4(j,:);
            b = alloc_vec_c4(i,:);
            diff = a - b;
            if ~all(diff >= 0)
                tran_matrix_c4(i, j)=0;
            else
                if b(1) == 0
                    tran_matrix_c4(i, j)=prob_LI(diff(1)+1);
                else
                    tran_matrix_c4(i, j)=1;
                end
            end
        end
    end

    tran_matrix{1,k}=tran_matrix_c1;
    tran_matrix{2,k}=tran_matrix_c2;
    tran_matrix{3,k}=tran_matrix_c3;
    tran_matrix{4,k}=tran_matrix_c4;
    
end


alloc_vec=cell(4,1);
alloc_vec{1,1}=alloc_vec_c1;
alloc_vec{2,1}=alloc_vec_c2;
alloc_vec{3,1}=alloc_vec_c3;
alloc_vec{4,1}=alloc_vec_c4;

%% Find index for state variable of each observation

state_index=zeros(N,1);

for i=1:N    
    temp_index=find(organ_available(i,:)~=0);    
    alloc_vec_temp=alloc_vec{category(i),1};    
    state_index(i)=find(ismember(alloc_vec_temp, organ_allocated(i, temp_index), 'rows'));    
end


%% Start MLE 

para_init=[10;5;3;0;0;0;0;0]; % initial guess of the parameters
options=optimset('Display','iter','MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-6,'TolX',1e-4);

gs = GlobalSearch;
llmin = @(para) log_ll(para, wait_cost, beta, T, N, ...
                    alloc_vec, category, tran_matrix, time_index, state_index, cross_clamp);
problem = createOptimProblem('fmincon', 'objective', llmin, 'x0', para_init);
[para_est, val] = run(gs, problem)

% para_est=fminsearch(@(para) log_ll(para, wait_cost, beta, T, N, ...
%                     alloc_vec, category, tran_matrix, time_index, state_index, cross_clamp), para_init, options)

save('para_est.mat','para_est','para_init');

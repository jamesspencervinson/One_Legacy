%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal Stopping Problem
% 05/03/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear

path='/Users/yixin/Dropbox/OneLegacy-Yi/matlab_code';
cd(path);


load('est_data_1h_v2.mat');
% time interval between two obs: 1 hour
% ignore intestine to reduce state space
% when cleaning the dataset: one donor has LU alone, dropped


%% load data and set parameters

organ_available=[HR_available, LU_available, LI_available];   
organ_allocated=[HR_allocated, LU_allocated, LI_allocated];


wait_cost=-1;
beta=1;
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

% HR
temp=zeros(N,2);
for i=1:N-1    
    temp(i, 1)=(organ_available(i, 1)==1 & organ_allocated(i, 1)==0 & donor_index(i)==donor_index(i+1));
    temp(i, 2)=(organ_available(i, 1)==1 & organ_allocated(i, 1)==0 & organ_allocated(i+1, 1)==1 & donor_index(i)==donor_index(i+1));
end
temp=sum(temp, 1);
prob_HR=[1-temp(2)/temp(1); 0; temp(2)/temp(1); 1];


% LU
temp=zeros(N,5);
for i=1:N-1    
    temp(i, 1)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==0 & donor_index(i)==donor_index(i+1));
    temp(i, 2)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==0 & organ_allocated(i+1, 2)==1 & donor_index(i)==donor_index(i+1));
    temp(i, 3)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==0 & organ_allocated(i+1, 2)==2 & donor_index(i)==donor_index(i+1));
    
    temp(i, 4)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==1 & donor_index(i)==donor_index(i+1));
    temp(i, 5)=(organ_available(i, 2)==2 & organ_allocated(i, 2)==1 & organ_allocated(i+1, 2)==2 & donor_index(i)==donor_index(i+1));
    
end
temp=sum(temp, 1);
prob_LU=[1-(temp(2)+temp(3))/temp(1); 0; 0; temp(2)/temp(1); 1-temp(5)/temp(4); 0; temp(3)/temp(1); temp(5)/temp(4); 1];



% LI
temp=zeros(N,2);
for i=1:N-1    
    temp(i, 1)=(organ_available(i, 3)==1 & organ_allocated(i, 3)==0 & donor_index(i)==donor_index(i+1));
    temp(i, 2)=(organ_available(i, 3)==1 & organ_allocated(i, 3)==0 & organ_allocated(i+1, 3)==1 & donor_index(i)==donor_index(i+1));
end
temp=sum(temp, 1);
prob_LI=[1-temp(2)/temp(1); 0; temp(2)/temp(1); 1];



%% Construct transition matrix


tran_matrix_c1=zeros(size(alloc_vec_c1,1), size(alloc_vec_c1,1));

for i=1:size(alloc_vec_c1,1)
    for j=1:size(alloc_vec_c1,1)
        
        temp1=prob_HR(ismember(HR_vec, [alloc_vec_c1(i, 1), alloc_vec_c1(j, 1)],'rows'));
        temp2=prob_LU(ismember(LU_vec, [alloc_vec_c1(i, 2), alloc_vec_c1(j, 2)],'rows'));
        temp3=prob_LI(ismember(LI_vec, [alloc_vec_c1(i, 3), alloc_vec_c1(j, 3)],'rows'));
        tran_matrix_c1(i, j)=temp1*temp2*temp3;
    end
end


tran_matrix_c2=zeros(size(alloc_vec_c2,1), size(alloc_vec_c2,1));

for i=1:size(alloc_vec_c2,1)
    for j=1:size(alloc_vec_c2,1)
        
        temp1=prob_HR(ismember(HR_vec, [alloc_vec_c2(i, 1), alloc_vec_c2(j, 1)],'rows'));
        temp2=prob_LI(ismember(LI_vec, [alloc_vec_c2(i, 2), alloc_vec_c2(j, 2)],'rows'));
        tran_matrix_c2(i, j)=temp1*temp2;
    end
end


tran_matrix_c3=zeros(size(alloc_vec_c3,1), size(alloc_vec_c3,1));

for i=1:size(alloc_vec_c3,1)
    for j=1:size(alloc_vec_c3,1)
        
        temp1=prob_LU(ismember(LU_vec, [alloc_vec_c3(i, 1), alloc_vec_c3(j, 1)],'rows'));
        temp2=prob_LI(ismember(LI_vec, [alloc_vec_c3(i, 2), alloc_vec_c3(j, 2)],'rows'));
        tran_matrix_c3(i, j)=temp1*temp2;
    end
end


tran_matrix_c4=zeros(size(alloc_vec_c4,1), size(alloc_vec_c4,1));

for i=1:size(alloc_vec_c4,1)
    for j=1:size(alloc_vec_c4,1)
        
        temp1=prob_LI(ismember(LI_vec, [alloc_vec_c4(i, 1), alloc_vec_c4(j, 1)],'rows'));
        tran_matrix_c4(i, j)=temp1;
    end
end


alloc_vec=cell(4,1);
alloc_vec{1,1}=alloc_vec_c1;
alloc_vec{2,1}=alloc_vec_c2;
alloc_vec{3,1}=alloc_vec_c3;
alloc_vec{4,1}=alloc_vec_c4;


tran_matrix=cell(4,1);
tran_matrix{1,1}=tran_matrix_c1;
tran_matrix{2,1}=tran_matrix_c2;
tran_matrix{3,1}=tran_matrix_c3;
tran_matrix{4,1}=tran_matrix_c4;


%% Find index for state variable of each observation

state_index=zeros(N,1);

for i=1:N    
    temp_index=find(organ_available(i,:)~=0);    
    alloc_vec_temp=alloc_vec{category(i),1};    
    state_index(i)=find(ismember(alloc_vec_temp, organ_allocated(i, temp_index), 'rows'));    
end


%% Start MLE 

para_init=[3;1;0.8]; % initial guess of the parameters
options=optimset('Display','iter','MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-4,'TolX',1e-4);


para_est=fminsearch(@(para) log_ll(para, wait_cost, beta, T, N, ...
                    alloc_vec, category, tran_matrix, time_index, state_index, cross_clamp), para_init, options); 

save('para_est.mat','para_init');

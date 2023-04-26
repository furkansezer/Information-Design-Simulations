%%random vector generation
n=5; %number of agents
h=10; %number of scenarios
pd = makedist('Normal');

rng('default') % For reproducibility
rand_numbers = random(pd,n,n,h); %

vec_norm=zeros(h,1);

for i=1:h
vec_norm(i)=norm(rand_numbers(:,:,i));
end

rand_vector=zeros(n,n,h);
for j=1:n
    for k=1:n
        for i=1:h
            rand_vector(j,k,i)=rand_numbers(j,k,i)/vec_norm(i);
        end  
    end    
end    

%%

eps1=0.1:0.1:1;
eps2=0.1;

% rho=1;


solution=zeros(h,length(eps1),2*n,2*n);

optimistic_solution=zeros(length(eps1),2*n,2*n);
pessimistic_solution=zeros(length(eps1),2*n,2*n);

opt_value=zeros(h,length(eps1));

max_opt_value=0;
min_opt_value=100;

S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

shift=zeros(2*n,2*n);

cov_theta=0.5*ones(n,n);

for i=1:n
    cov_theta(i,i)=5;
end 

    
Q=-1*ones(n);
for i=1:n
    Q(i,i)=5;
end 

no_info_solution=[zeros(n,n),zeros(n,n);zeros(n,n),Q];
full_info_solution=[Q\cov_theta/inv(Q).',Q\cov_theta;cov_theta/inv(Q).',cov_theta];

V_ss=zeros(2*n,2*n); %agreement
    for i=1:n
        for j=1:n
            V_ss(i,j)= 1/n;
        end 
            V_ss(i,i)=-1 + 1/n;
    end

V_sw=[-Q eye(n);eye(n) zeros(n)];    %social welfare
    
for g=1:h % g for scenarios



for w=1:10 % w for epsilon variation

overall_shift=[eps2*ones(n) zeros(n);zeros(n) zeros(n)];
for i=1:n
        overall_shift(i,i)=eps1(w);
end   


perturbed_h=overall_shift(1:5,1:5).*squeeze(rand_vector(:,:,g));

perturbed_big=[perturbed_h zeros(n);zeros(n) zeros(n)];
    
cvx_begin quiet
variable X(2*n,2*n) symmetric semidefinite
variable t

maximize (t)
subject to

% trace((V_ss*X.'))>=t; %aggreement 

trace((V_sw*X.'))-trace(perturbed_big*X')>=t; %social welfare


for i=1:n
    for j=1:i
        if(i==j)
            S(n+i,n+j)=1;
            S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
            S(n+i,n+j)=0;
            
        end
        if (j<i)
            S(n+i,n+j)=1;
           S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
       
             S(n+i,n+j)=0;
        end
    end
end

for k=1:n
    for i=1:n
        R(k,i)=Q(k,i)/2;
        R(i,k)=Q(k,i)/2;
        shift(k,i)=rand_vector(k,i,g)*(eps2/2);
        shift(i,k)=rand_vector(i,k,g)*(eps2/2);
                
    end
    R(k,k)=Q(k,k);
    shift(k,k)=eps1(w)*rand_vector(k,k,g);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    
    
    trace(shift*X')+trace(R*X')==0; % second constraint group   
    

    R=zeros(2*n,2*n);
    shift=zeros(2*n,2*n);  
    
end


% X>=0; %third constraint


cvx_end 

solution(g,w,:,:)=X;
opt_value(g,w)=cvx_optval;

end

% if cvx_optval>max_opt_value
%     max_opt_value=cvx_optval;
%     optimistic_solution=X;
% end 
% 
% if cvx_optval<min_opt_value
%     min_opt_value=cvx_optval;
%     pessimistic_solution=X;
% end 

end

 for g=1:10
     plot(eps1, opt_value(g,:), 'LineWidth',2)
     hold on
 end    
 hold off

    
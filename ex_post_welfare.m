a=5; %prior parameters
b=1;

n=30; %number of initial random gamma priors considered
n1=10; %number of datapoints used to update gamma prior
n2=10; %number of samples from each gamma posterior
true_mean=1/4;
w_s=1/5;

c=0.2;%marginal cost


suff_full=zeros(n,1);
suff_no=zeros(n,1);
posterior_full_sample=zeros(n,n1);
posterior_no_sample=zeros(n,n1);

prior=makedist('Gamma','a',3,'b',1); 
prior_sample=random(prior,n,1);

likelihood=makedist('Exponential','mu',true_mean); %real likelihood
likelihood_seller=makedist('Exponential','mu',w_s); %real likelihood

for i=1:n
    
suff_full(i)=sum(random(likelihood,n1,1)); %sufficient statistic under full information
suff_no(i)= prior_sample(i)*n1; %sufficient statistic under no information

posterior_full=makedist('Gamma','a',a+n1,'b',b+suff_full(i));
posterior_no=makedist('Gamma','a',a+n1,'b',b+suff_no(i));

posterior_full_sample(i,:)=random(posterior_full,1,n2);
posterior_no_sample(i,:)=random(posterior_no,1,n2);
end

%Ex-post prices
realized_benefit=random(likelihood,n,n2); %b^{B}
realized_benefit_seller=random(likelihood_seller,n,n2);  %b^{S}
price_full=zeros(n,n2);
price_no=zeros(n,n2);
price_seller=zeros(n,n2);

for i=1:n
    for j=1:n2
price_full(i,j)=((posterior_full_sample(i,j)*c+2)*exp(-posterior_full_sample(i,j)*realized_benefit(i,j))+c+2)/(2*posterior_full_sample(i,j)*exp(-posterior_full_sample(i,j)*realized_benefit(i,j)));
price_no(i,j)=((posterior_no_sample(i,j)*c+2)*exp(-posterior_no_sample(i,j)*realized_benefit(i,j))+c+2)/(2*posterior_no_sample(i,j)*exp(-posterior_no_sample(i,j)*realized_benefit(i,j)));
price_seller(i,j)=((w_s*c+2)*exp(-w_s*realized_benefit_seller(i,j))+c+2)/(2*w_s*exp(-w_s*realized_benefit_seller(i,j)));
    
    end
end

%Welfare calculation
welfare_full=zeros(n,n2);
welfare_no=zeros(n,n2);

for i=1:n
    for j=1:n2
welfare_full(i,j)=exp(-posterior_full_sample(i,j)*realized_benefit(i,j))*exp(-w_s*price_seller(i,j))/w_s+exp(-w_s*realized_benefit_seller(i,j))*exp(-posterior_full_sample(i,j)*price_full(i,j))/posterior_full_sample(i,j);
welfare_no(i,j)=exp(-posterior_no_sample(i,j)*realized_benefit(i,j))*exp(-w_s*price_seller(i,j))/w_s+exp(-w_s*realized_benefit_seller(i,j))*exp(-posterior_no_sample(i,j)*price_no(i,j))/posterior_no_sample(i,j);

    end
end

%Welfare comparison
welfare_delta=welfare_full-welfare_no;
mean_welfare_delta=mean(welfare_delta,2);
mean_full_welfare=mean(welfare_full,2);
mean_no_welfare=mean(welfare_no,2);

mean_price_full_buyers=mean(price_full,2);
mean_price_no_buyers=mean(price_no,2);

mean_price_sellers=mean(price_seller,2);
scatter(log(mean_price_full_buyers),mean_welfare_delta,100,'filled') %recommended prices under full info and welfare difference
scatter(log(mean_price_full_buyers),mean_full_welfare)
scatter(log(mean_price_no_buyers),mean_no_welfare)

mean_benefits_buyers=mean(realized_benefit,2);
mean_benefits_sellers=mean(realized_benefit_seller,2);

% 
% scatter(mean_benefits_buyers,mean_full_welfare)
% scatter(mean_benefits_buyers,mean_no_welfare)
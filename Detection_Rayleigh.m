clc;
clear all;

% APs Jointly estimate

MRC_emperical = [];
MRC_theoretical = [];
MRC_emperical_equation = [];

user_location_total = [];

for t = 1:1000
    
    cell_size = 500;
    shadow = 8;
    M = 50;
    N = 2;
    K = 100;
    K_active = 10;
    tau = 5;
    SNR_db = 111;
    SNR_linear = 10^(SNR_db/10);
    M = 30:5:50;
    
    parfor i = 1:size(M,2)
        
        bit = 2*ones(M(i),1);
        alpha = 1-pi*sqrt(3)/2*2.^(-2*bit);
        
        success_rate_CF_K1(t,i) = ADMM(cell_size,shadow,M(i),N,K,10,tau,SNR_linear,alpha,'CF');
        success_rate_CF_K2(t,i) = ADMM(cell_size,shadow,M(i),N,K,20,tau,SNR_linear,alpha,'CF');
        
        bit = 2*ones(1,1);
        alpha = 1-pi*sqrt(3)/2*2.^(-2*bit);
        
        success_rate_CO_K1(t,i) = ADMM(cell_size,shadow,1,M(i)*N,K,10,tau,SNR_linear,alpha,'CO');
        success_rate_CO_K2(t,i) = ADMM(cell_size,shadow,1,M(i)*N,K,20,tau,SNR_linear,alpha,'CO');
    end
    
end

plot(M,mean(success_rate_CF_K1),'LineWidth',2,'Color','#D95319','LineStyle','--'); hold on;
plot(M,mean(success_rate_CF_K2),'LineWidth',2,'Color','#A2142F','LineStyle','-');
plot(M,mean(success_rate_CO_K1),'LineWidth',2,'Color','#4DBEEE','LineStyle',':');
plot(M,mean(success_rate_CO_K2),'LineWidth',2,'Color','#0072BD','LineStyle','-.');
legend('distributed ($K^a=10$)','distributed ($K^a=20$)','co-located ($K^a=10$)','co-located ($K^a=20$)','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
xlabel('$M$','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Detection rate','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');

function [success_rate] = ADMM(cell_size,shadow,M,N,K,K_active,tau,SNR_linear,alpha,key)
% Generate channel
h = cellfun(@(~)(randn(N,K)+1i*randn(N,K))/sqrt(2),cell(M,1),'UniformOutput',false);
if strcmp(key,'CF')
    base_location = (rand(M,1)-0.5)*2*cell_size+1i*(rand(M,1)-0.5)*2*cell_size;
elseif strcmp(key,'CO')
    base_location = zeros(M,1);
end
user_location = (rand(K,1)-0.5)*2*cell_size+1i*(rand(K,1)-0.5)*2*cell_size;
beta = (abs(base_location-user_location.')/1).^(-3.5).*10.^(10^(shadow/10)*(randn(M,K))/10);
g = cellfun(@(a,b)repmat(a,N,1).*b,mat2cell(sqrt(beta),ones(M,1),K),h,'UniformOutput',false);

% Transmit pilot
pilot = (normc(randn(tau,K))+1i*normc(randn(tau,K)))/sqrt(2);
x_p = pilot;
x_p = x_p./repmat(sum(abs(x_p).^2),tau,1);
% Activity
index = randperm(K,K_active); activity = zeros(K,1); activity(index) = 1;
y_p_clean=cellfun(@(a)sqrt(SNR_linear*tau)*a*diag(activity)*x_p',g,'UniformOutput',false);
y_p = cellfun(@(a)a+(randn(N,tau)+1i*randn(N,tau))/sqrt(2),y_p_clean,'UniformOutput',false);

% Transmit pilot after quantization
y_p = cellfun(@(a,b,c)a*b+(mvnrnd(zeros(N,1),c*eye(N),tau)+1i*mvnrnd(zeros(N,1),c*eye(N),tau))'/sqrt(2),...
    y_p,num2cell(alpha,M),num2cell(alpha.*(1-alpha).*(SNR_linear*sum(beta.*repmat(activity',M,1),2)+1)),'UniformOutput',false);

estimated_activity = mat2cell(zeros(K,M),K,ones(M,1))';
multiplier = mat2cell(zeros(K,M),K,ones(M,1))';
common_estimated_activity = zeros(K,1);

for outer_iteration = 1:10
    estimated_activity = cellfun(@(a,b,c,d,e)update_step(N,K,tau,SNR_linear,a,b,activity,x_p,c,d,common_estimated_activity,e),...
        num2cell(alpha),mat2cell(beta,ones(M,1),K),y_p,estimated_activity,multiplier,'UniformOutput',false);
    common_estimated_activity = mean(cell2mat(estimated_activity')+cell2mat(multiplier'),2);
    multiplier = cellfun(@(a,b)a+b-common_estimated_activity,multiplier,estimated_activity,'UniformOutput',false);
    temp(outer_iteration) = max(cell2mat(estimated_activity')-common_estimated_activity,[],'all');
    if temp(outer_iteration)<1e-2
        break;
    end
end

[~,found_index]=maxk(common_estimated_activity,K_active);
success_rate = 1-size(setdiff(index,found_index),2)/K_active;
end

function [estimated_activity] = update_step(N,K,tau,SNR_linear,alpha,beta,activity,x_p,y_p,estimated_activity,common_estimated_activity,multiplier)
for inner_iteration = 1:10
    estimated_activity_old = estimated_activity;
    for k = 1:K
        estimated_activity_temp = estimated_activity; estimated_activity_temp(k) = 0;
        sigma_noise = alpha^2 + alpha*(1-alpha)*(SNR_linear*sum(beta.*activity')+1);
        % Necessary values
        covariance_k = inv(tau*SNR_linear*alpha^2*x_p*diag(estimated_activity_temp.*beta')*x_p' + sigma_noise*eye(tau));
        s = real(x_p(:,k)'*(covariance_k)*x_p(:,k));
        q = x_p(:,k)'*(covariance_k)*y_p'/sqrt(N);
        % ADMM
        z = common_estimated_activity(k);
        v = multiplier(k);
        p = tau*SNR_linear*alpha^2*beta(k);
        potential=roots([p^2*s^2, ((v-z)*p^2*s^2+2*p*s), (p^2*s^2+2*p*s*(v-z)+1), -p*(q*q')+v-z+p*s]);
        estimated_activity(k)=max([potential(imag(potential(:))==0)',0]);
    end
    temp(inner_iteration) = max(abs(estimated_activity-estimated_activity_old));
    if temp(inner_iteration) < 1e-2
        break;
    end
end
end
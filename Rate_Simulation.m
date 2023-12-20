clc;
clear all;

MRC_emperical_CF = [];
MRC_theoretical_CF = [];

user_location_total = [];

for t = 1:1000
    
    size = 500;
    shadow = 8;
    M = 50;
    N = 2;
    K = 10;
    tau = 5;
    SNR_db = 111;
    SNR_linear = 10^(SNR_db/10);
    bit = 2*ones(M,1);
    alpha = 1-pi*sqrt(3)/2*2.^(-2*bit);

    [MRC_emperical_CF(:,t),MRC_theoretical_CF(:,t)] = Rate(size,shadow,M,N,K,tau,SNR_linear,alpha,'CF');

    % Colocated
    N = M*N; M = 1; bit = 2*ones(M,1); alpha = 1-pi*sqrt(3)/2*2.^(-2*bit);
    [MRC_emperical_CO(:,t),MRC_theoretical_CO(:,t)] = Rate(size,shadow,M,N,K,tau,SNR_linear,alpha,'CO'); 
end

figure(1);
h1=cdfplot(MRC_theoretical_CO(:)); hold on;
h2=cdfplot(MRC_emperical_CO(:)); hold on;
h3=cdfplot(MRC_theoretical_CF(:)); hold on;
h4=cdfplot(MRC_emperical_CF(:)); hold on;
legend('co-located (analytical)','co-located (numerical)','distributed (analytical)','distributed (numerical)','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
axis([0 16 0 1]);
xlabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('CDF','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
title('')
set(h1,'LineWidth',2,'Color','#0072BD','LineStyle','-.');
set(h2,'LineWidth',2,'Color','#4DBEEE','LineStyle',':');
set(h3,'LineWidth',2,'Color','#A2142F','LineStyle','-');
set(h4,'LineWidth',2,'Color','#D95319','LineStyle','--');

figure(2)
subplot(2,1,1)
histogram(MRC_emperical_CO(:),'Normalization','probability')
axis([0 16 0 1.0]);
xlabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('PMF','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
annotation('textbox',[0.66 0.80 .18 .08],'String','distributed','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
subplot(2,1,2)
histogram(MRC_emperical_CF(:),'Normalization','probability')
axis([0 16 0 0.1]);
xlabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('PMF','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
annotation('textbox',[0.66 0.33 .18 .08],'String','co-located','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman');

function [MRC_emperical,MRC_theoretical] = Rate(size,shadow,M,N,K,tau,SNR_linear,alpha,key)

% Generate channel
h = cellfun(@(~)(randn(N,K)+1i*randn(N,K))/sqrt(2),cell(M,1),'UniformOutput',false);
if strcmp(key,'CF')
    base_location = (rand(M,1)-0.5)*2*size+1i*(rand(M,1)-0.5)*2*size;
elseif strcmp(key,'CO')
    base_location = zeros(M,1);
end
user_location = (rand(K,1)-0.5)*2*size+1i*(rand(K,1)-0.5)*2*size;
beta = (abs(base_location-user_location.')/1).^(-3.5).*10.^(10^(shadow/10)*(randn(M,K))/10);
g = cellfun(@(a,b)repmat(a,N,1).*b,mat2cell(sqrt(beta),ones(M,1),K),h,'UniformOutput',false);

% Transmit pilot
x_p = (normc(randn(tau,K))+1i*normc(randn(tau,K)))/sqrt(2);
x_p = x_p./repmat(sum(abs(x_p).^2),tau,1);
y_p_clean=cellfun(@(a)sqrt(SNR_linear*tau)*a*x_p',g,'UniformOutput',false);
y_p = cellfun(@(a)a+(randn(N,tau)+1i*randn(N,tau))/sqrt(2),y_p_clean,'UniformOutput',false);

% Transmit pilot after quantization
y_p = cellfun(@(a,b,c)a*b+(mvnrnd(zeros(N,1),c*eye(N),tau)+1i*mvnrnd(zeros(N,1),c*eye(N),tau))'/sqrt(2),...
    y_p,num2cell(alpha,M),num2cell(alpha.*(1-alpha).*(SNR_linear*sum(beta,2)+1)),'UniformOutput',false);

% MMSE channel estimation
g_estimation = cellfun(@(a,b,c)a*x_p*diag((sqrt(SNR_linear*tau)*b)./(SNR_linear*((1-c)*sum(b)+c*tau*b*abs(x_p'*x_p).^2)+1)),y_p,mat2cell(beta,ones(M,1),K),num2cell(alpha),'UniformOutput',false);
gamma = cellfun(@(b,c)((c*SNR_linear*tau*b.^2)./(SNR_linear*((1-c)*sum(b)+c*tau*b*abs(x_p'*x_p).^2)+1)),mat2cell(beta,ones(M,1),K),num2cell(alpha),'UniformOutput',false);

% Transmit data
x_d = (randn(K,1)+1i*randn(K,1))/sqrt(2); x_d = x_d./abs(x_d);
y_d_clean=cellfun(@(a)sqrt(SNR_linear)*a*x_d,g,'UniformOutput',false);
y_d = cellfun(@(a)a+(randn(N,1)+1i*randn(N,1))/sqrt(2),y_d_clean,'UniformOutput',false);

% Transmit data after quantization
y_d = cellfun(@(a,b,c)a*b+(mvnrnd(zeros(N,1),c*eye(N),1)+1i*mvnrnd(zeros(N,1),c*eye(N),1))'/sqrt(2),...
    y_d,num2cell(alpha,M),num2cell(alpha.*(1-alpha).*(SNR_linear*sum(beta,2)+1)),'UniformOutput',false);

% MRC detection
% Emperical
MRC_emperical=(log2(1+abs(x_d).^2./abs(sum(cell2mat(cellfun(@(a,b)a'*b,g_estimation,y_d,'UniformOutput',false)'),2)./(sqrt(SNR_linear)*N*(alpha'*cell2mat(gamma))')-x_d).^2));

% Theoretical
% First part in denominator
for k = 1:K
    tempi = 0;
    for i = 1:K
        if i~=k
            tempm=0;
            for m = 1:M
                tempm = tempm + alpha(m)*gamma{m}(k)/beta(m,k)*beta(m,i);
            end
            tempi = tempi + tempm^2*abs(x_p(:,k)'*x_p(:,i))^2;
        end
    end
    first(k,1) = SNR_linear*N^2*tempi;
end

second = N*SNR_linear*sum(cell2mat(gamma)'*(beta.*repmat(alpha.^2,1,K)),2);
% Third part in denominator
third = N*((alpha.*(1+(1-alpha)*SNR_linear.*sum(beta,2)))'*cell2mat(gamma))';
fourth = (SNR_linear*N^2*((ones(M,1))'*cell2mat(gamma)).^2)'-(SNR_linear*N^2*((alpha)'*cell2mat(gamma)).^2)';

MRC_theoretical=(log2(1+(SNR_linear*N^2*(alpha'*cell2mat(gamma)).^2)'./(first+second+third)));
end
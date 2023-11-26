clc;
clear all;



% MRC_emperical_total_error = [];
% MRC_theoretical_total_error = [];
% MMSE_emperical_total_error = [];
% MMSE_theoretical_total_error = [];

MRC_emperical_CF = [];
MRC_theoretical_CF = [];
% MRC_emperical_equation = [];
% Imperfect_MMSE_theoretical = [];

user_location_total = [];

for t = 1:1000

size = 500;
shadow = 8;
M = 100;
N = 20;
K = 20;
tau = 50;
SNR_db = 111;
SNR_linear = 10^(SNR_db/10);
bit = 2*ones(M,1);
alpha = 1-pi*sqrt(3)/2*2.^(-2*bit);%%%%%ones(M,1);%



% figure(1);
% plot(base_location,'ro');plot(user_location,'bd');
[MRC_emperical_CF(:,t),MRC_theoretical_CF(:,t)] = Rate(size,shadow,M,N,K,tau,SNR_linear,alpha,'CF');



% Colocated
[MRC_emperical_CO(:,t),MRC_theoretical_CO(:,t)] = Rate(size,shadow,M,N,K,tau,SNR_linear,alpha,'CO');





end


% signal = 0.1;
% noise=20*1e6 * 1.381*1e-23 * 290 * 10^(9/10);
% 10*log10(signal/noise)
% signal/noise

% mean(MRC_emperical_total_error,'all')
% mean(MRC_theoretical_total_error,'all')
% cdfplot(MRC_emperical_total_error(:)); hold on;
% cdfplot(MRC_theoretical_total_error(:))
% legend('MRC emperical','MRC theoretical','MMSE emperical','MMSE theoretical')

% mean(MMSE_emperical_total_error,'all')
% mean(MMSE_theoretical_total_error,'all')
% cdfplot(MMSE_emperical_total_error(:)); hold on;
% cdfplot(MMSE_theoretical_total_error(:))
% legend('MMSE emperical','MMSE theoretical')


% 10*log10((100/(20*10^6*1.381*10^(-23)*290*10^(9/10)*1000)))

figure(1);
% mean(MRC_emperical,'all')
% mean(MRC_emperical_equation,'all')
% mean(MRC_theoretical,'all')
h1=cdfplot(MRC_emperical_CF(:)); hold on;
h2=cdfplot(MRC_theoretical_CF(:)); hold on;
h3=cdfplot(MRC_emperical_CO(:)); hold on;
h4=cdfplot(MRC_theoretical_CO(:)); hold on;
legend('distributed (simulation)','distributed (theory)','co-located (simulation)','co-located (theory)','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
axis([0 20 0 1]);
xlabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('CDF','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
title('')
set(h1,'LineWidth',2,'Color','#D95319','LineStyle','--');
set(h2,'LineWidth',2,'Color','#A2142F','LineStyle','-');
set(h3,'LineWidth',2,'Color','#4DBEEE','LineStyle',':');
set(h4,'LineWidth',2,'Color','#0072BD','LineStyle','-.');



% mean(Imperfect_MMSE_emperical,'all')
% mean(Imperfect_MMSE_theoretical,'all')
% h3=cdfplot(Imperfect_MMSE_emperical(:)); hold on;
% h4=cdfplot(Imperfect_MMSE_theoretical(:)); hold on;
% legend('MRC emperical','MRC theoretical','MMSE emperical', 'MMSE theoretical','interpreter','latex','FontSize', 18, 'FontName', 'Times New Roman')




% set(h3,'LineWidth',2)
% set(h4,'LineWidth',2)
% 
figure(3)
subplot(2,1,1)
histogram(MRC_emperical_CO(:),'Normalization','probability')
axis([0 20 0 0.6]);
xlabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('PMF','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
annotation('textbox',[0.66 0.80 .18 .08],'String','distributed','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
subplot(2,1,2)
histogram(MRC_emperical_CF(:),'Normalization','probability')
axis([0 20 0 0.06]);
xlabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('PMF','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
annotation('textbox',[0.66 0.33 .18 .08],'String','co-located','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')




% figure(4)
% [B,I] = mink(Imperfect_MMSE_emperical,0.01*K*1000);
% plot(user_location_total(I),'bd')
% axis([-500 500 -500 500]);

function [MRC_emperical_error,MRC_theoretical_error] = MRC(x_p,y_p,tau,SNR_linear,beta,g)
MRC_emperical_error=mean(abs(x_p'*y_p/sqrt(SNR_linear*tau)-g.').^2,'all')/mean(abs(g).^2,'all');
MRC_theoretical_error=mean((beta*abs(x_p'*x_p-diag(diag(x_p'*x_p))).^2+1/(SNR_linear*tau)),'all')/mean(abs(g).^2,'all');
end

function [MMSE_emperical_error,MMSE_theoretical_error] = MMSE(x_p,y_p,tau,SNR_linear,beta,g)
MMSE_emperical_error=sum(abs(sqrt(SNR_linear*tau)*beta.*((x_p'*y_p).')./(SNR_linear*tau*beta*abs(x_p'*x_p).^2+1)-g).^2./repmat(sum(abs(g).^2),size(g,1),1));
MMSE_theoretical_error=sum((beta-SNR_linear*tau*beta.^2./(SNR_linear*tau*beta*abs(x_p'*x_p).^2+1))./repmat(sum(abs(g).^2),size(g,1),1));
end

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
ideal_pilot = diag(ones(tau,1)+1i*ones(tau,1))/sqrt(2); index = randperm(tau,K);ideal_pilot(:,index);%
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



% g_estimation = cellfun(@(a,b)a*x_p*diag((sqrt(SNR_linear*tau)*b)./(SNR_linear*tau*b*abs(x_p'*x_p).^2+1)),y_p,mat2cell(beta,ones(M,1),K),'UniformOutput',false);
% gamma=cellfun(@(a,b)(a*SNR_linear*tau*b.^2)./(SNR_linear*tau*b*abs(x_p'*x_p).^2+1),num2cell(alpha),mat2cell(beta,ones(M,1),K),'UniformOutput',false);



% t1=cell2mat(cellfun(@(a,b)mean(abs(a-b).^2),g_estimation,g,'UniformOutput',false));
% t2=cell2mat(cellfun(@(a,b)(a-b),mat2cell(beta,ones(M,1),K),gamma,'UniformOutput',false));
% subplot(2,1,1);histogram(t1(:),200);subplot(2,1,2);histogram(t2(:),200);
% 
% mean(t1)
% mean(t2)


% Transmit data
x_d = (randn(K,1)+1i*randn(K,1))/sqrt(2); x_d = x_d./abs(x_d);
y_d_clean=cellfun(@(a)sqrt(SNR_linear)*a*x_d,g,'UniformOutput',false);
y_d = cellfun(@(a)a+(randn(N,1)+1i*randn(N,1))/sqrt(2),y_d_clean,'UniformOutput',false);

y_d_original = y_d;

% Transmit data after quantization
y_d = cellfun(@(a,b,c)a*b+(mvnrnd(zeros(N,1),c*eye(N),1)+1i*mvnrnd(zeros(N,1),c*eye(N),1))'/sqrt(2),...
    y_d,num2cell(alpha,M),num2cell(alpha.*(1-alpha).*(SNR_linear*sum(beta,2)+1)),'UniformOutput',false);


% MRC detection
% Emperical
MRC_emperical=(log2(1+abs(x_d).^2./abs(sum(cell2mat(cellfun(@(a,b)a'*b,g_estimation,y_d,'UniformOutput',false)'),2)./(sqrt(SNR_linear)*N*(alpha'*cell2mat(gamma))')-x_d).^2));


% % Emperical equation
% A = (SNR_linear*N^2*(alpha'*cell2mat(gamma)).^2)';
% B = abs(sum(cell2mat(cellfun(@(a,b,c)sqrt(SNR_linear)*c*diag(a'*b),g_estimation,g,num2cell(alpha),'UniformOutput',false)'),2)-sqrt(A)).^2;
% C_temp = cellfun(@(a,b,c)sqrt(SNR_linear)*c*(a'*b-diag(diag(a'*b))),g_estimation,g,num2cell(alpha),'UniformOutput',false);
% C=sum(abs(sum(cat(3,C_temp{:}),3)).^2,2);
% D = abs(sum(cell2mat(cellfun(@(a,b,c,d)a'*(b-c)*d,g_estimation,y_d_original, y_d_clean,num2cell(alpha),'UniformOutput',false)'),2)).^2;
% E = abs(sum(cell2mat(cellfun(@(a,b,c,d)a'*(c-b*d),g_estimation,y_d_original,y_d,num2cell(alpha),'UniformOutput',false)'),2)).^2;
% MRC_emperical_equation=[MRC_emperical_equation;log2(1+A./(B+C+D+E))];


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
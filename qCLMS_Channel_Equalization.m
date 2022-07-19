%% Channel Equalization using Least Mean Square (LMS) algorithm
% Author: SHUJAAT KHAN
% cite: Sadiq, Alishba, Imran Naseem, Shujaat Khan, Muhammad Moinuddin, Roberto Togneri, and Mohammed Bennamoun. "A novel quantum calculus-based complex least mean square algorithm (q-CLMS)." Applied Intelligence (2022): 1-20.

clc;
clear all;
close all;
pskmod_true = 1; % same as in paper
%pskmod_true = 0; % In case if you don't have Communication Toolbox

%% Channel and noise level

%h = [0.9 0.3 0.5 -0.1]; % Channel
SNRr = 20;              % Noise Level

%% LMS parameters
runs = 100;        % Number of independent runs (training repetation)
eta_CLMS = 4e-2;         % Learning rate / step size
eta_qCLMS = 1e-2;
order=20;           % Order of the equalizer

%% Visualization settings
f_size=14;   % Font size of figure
lw=2;       % linewidth of plot
f_FCLMS = [0.9 0.75];
mu_FCLMS = 1e-2;%9e-4;
mu_FoCLMS = 5e-2;
q=ones(order,1);
eU=zeros(order,1);
eta_max = 0.3;
%% Algorithm
for run = 1 : runs
    %% intialize weights
    h = 1 + 0*complex(randn(5,1),randn(5,1));
    U = zeros(order,1); % Input frame
    W = randn(order,1); % Initial Weigths
    Wq = W;
    W_FCLMS1 = W;
    W_FCLMS2 = W; 
    W_FoCLMS1 = W;
    W_FoCLMS2 = W;
    %% Input/Output data
    N = 300;               % Number of samples
    Bits = 16;               % Number of bits for modulation (2-bit for Binary modulation)
    if(pskmod_true)
        data = randi([0 1],1,N);        % Random signal
        d = real(pskmod(data,Bits));    % BPSK Modulated signal (desired/output)
    else
        %d = randn(1,N);                 % Gaussian Random signal
        d = randi([0 2^(Bits-1)],1,N)./(2^(Bits-1));        % Uniform Random signal
    end
    r = filter(h,1,d);              % Signal after passing through channel
    x = awgn(r, SNRr);              % Noisy Signal after channel (given/input)

    for n = 1 : N
        U(1,2:end) = U(1,1:end-1);  % Sliding window
        U(1,1) = x(n);              % Present Input
        %CLMS
        y = (W')*U;             % Calculating output of LMS
        e = d(n) - y;           % Instantaneous error 
        W = W +  0.05*eta_CLMS * conj(e) * U ;  % Weight update rule of LMS
        CLMS(run,n) = e * e';        % Instantaneous square error
        
        % Eq-CLMS
        Yq = (Wq')*U;
        eq = d(n)-Yq;
        eU(1,2:end) = eU(1,1:end-1);
        eU(1,1) = abs(eq);
        
        if((abs(eq))>(eta_max)) % original
        eU(:,:) = 1*eta_max;
        end
        
        q = (eU);
        
        Wq = Wq + 1.5*eta_qCLMS*conj(eq)*U.*q;
        EqCLMS(run,n) = eq * eq';
        
        %CFLMS
          y_FCLMS1 = (W_FCLMS1')*U;
        e_FCLMS1 = d(n) - y_FCLMS1;      

        W_FCLMS1 = W_FCLMS1 + 0.1*mu_FCLMS*conj(e_FCLMS1)*U ...,
            +  (0.1*mu_FCLMS/gamma(2-f_FCLMS(1)))*conj(e_FCLMS1)*U.*(W_FCLMS1.^(1-f_FCLMS(1)));
         K_FCLMS(run,n) = e_FCLMS1 * e_FCLMS1';
         
         
         y_FCLMS2 = (W_FCLMS2')*U;
        e_FCLMS2 = d(n) - y_FCLMS2;      

        W_FCLMS2 = W_FCLMS2 + 0.1*mu_FCLMS*conj(e_FCLMS2)*U ...,
            +  (0.1*mu_FCLMS/gamma(2-f_FCLMS(2)))*conj(e_FCLMS2)*U.*(W_FCLMS2.^(1-f_FCLMS(2)));
         K_FCLMS2(run,n) = e_FCLMS2 * e_FCLMS2';
         
         %FoCLMS
         
          y_FoCLMS1 = (W_FoCLMS1')*(U);
        e_FoCLMS1 = d(n) - y_FoCLMS1;      
A = 0.05*mu_FoCLMS/(2*gamma(2-f_FCLMS(1)));
      W_FoCLMS1 = W_FoCLMS1 + A*(conj(e_FoCLMS1)*U).*(real(W_FoCLMS1).^(1-f_FCLMS(1))+imag(W_FoCLMS1).^(1-f_FCLMS(1)))...,
        + A*(e_FoCLMS1*conj(U)).*(real(W_FoCLMS1).^(1-f_FCLMS(1))-imag(W_FoCLMS1).^(1-f_FCLMS(1)));
                
        Q_FoCLMS(run,n) = e_FoCLMS1* e_FoCLMS1';
        
          y_FoCLMS2 = (W_FoCLMS2')*(U);
        e_FoCLMS2 = d(n) - y_FoCLMS2;      
A = 0.05*mu_FoCLMS/(2*gamma(2-f_FCLMS(2)));
      W_FoCLMS2 = W_FoCLMS2 + A*(conj(e_FoCLMS2)*U).*(real(W_FoCLMS2).^(1-f_FCLMS(2))+imag(W_FoCLMS2).^(1-f_FCLMS(2)))...,
        + A*(e_FoCLMS2*conj(U)).*(real(W_FoCLMS2).^(1-f_FCLMS(2))-imag(W_FoCLMS2).^(1-f_FCLMS(2)));
                
        Q_FoCLMS2(run,n) = e_FoCLMS2* e_FoCLMS2';

        
        
    end
end

%% Calculation of performance parameters

MSE_CLMS = mean(CLMS,1);     % Mean square error
MSE_EqCLMS = mean(EqCLMS,1);
MSE_FCLMS = mean(K_FCLMS,1);
MSE_FCLMS2 = mean(K_FCLMS2,1);
MSE_FoCLMS = mean(Q_FoCLMS,1);
MSE_FoCLMS2 = mean(Q_FoCLMS2,1);

%% Plots
figure % MSE
% plot(10*log10(MJ),'->k','linewidth',lw)
plot(10*log10(MSE_CLMS),'b','linewidth',lw)
hold on
plot(10*log10(MSE_FCLMS),'cy','linewidth',lw)
plot(10*log10(MSE_FCLMS2),'cy--','linewidth',lw)
plot(10*log10(MSE_FoCLMS),'mag','linewidth',lw)
plot(10*log10(MSE_FoCLMS2),'mag--','linewidth',lw)
plot(10*log10(MSE_EqCLMS),'k','linewidth',lw)
lgd = legend('CLMS','CFLMS(0.9)','CFLMS(0.75)','FoCLMS(0.9)','FoCLMS(0.75)','Eq-CLMS')
xlabel('Number of iterations','FontSize',14,'FontWeight','bold','Color','k')
ylabel('MSE (dB)','FontSize',14,'FontWeight','bold','Color','k')
set(lgd,'FontSize',f_size)
grid minor
% save('Results_ChannelEqualization40.mat')


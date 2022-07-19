%% Author: Shujaat Khan
function [ y ] = awgnSK( x,snr_db )
%AWUN This function will add Gaussian distributed white noise of defined 
% SNR value of (snr_db) into the given signal (x).
SP=norm(x).^2; % Input signal power
NP=sqrt(SP*10^-(snr_db/10)/length(x)); % Desired noise level
y=x+NP*(randn(size(x))); % adding noise
end

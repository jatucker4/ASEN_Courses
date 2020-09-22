function P1 = one_sided_FFT(Y)
%ONE_SIDED_FFT Summary of this function goes here
%   Detailed explanation goes here
% NFFT= 2^(nextpow2(length(x))); 
% % Take fft, padding with zeros so that length(FFTX) is equal to NFFT 
% FFTX = fft(x,NFFT); 
% % Calculate the numberof unique points 
% NumUniquePts = ceil((NFFT+1)/2); 
% % FFT is symmetric, throw away second half 
% FFTX = FFTX(1:NumUniquePts); 
% % Take the magnitude of fft of x 
% MX = abs(FFTX); 
% % Scale the fft so that it is not a function of the length of x 
% % MX = MX/length(x); 
% % % Take the square of the magnitude of fft of x. 
% % MX = MX.^2; 
% % % Multiply by 2 because you threw out second half of FFTX above 
% % MX = MX*2; 
% % % DC Component should be unique. 
% % MX(1) = MX(1)/2; 
% % % Nyquist component should also be unique.
% % if ~rem(NFFT,2) 
% %    % Here NFFT is even; therefore,Nyquist point is included. 
% %    MX(end) = MX(end)/2; 
% % end 
L = length(Y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
end


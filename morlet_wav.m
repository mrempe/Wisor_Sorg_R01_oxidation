% ---- Chris Fink kindly provided this script in July 2019 ----
% -------------------------------------------------------------
% ---- The comments below are his -----------------------------
%
%this function performs time-frequency analysis on a signal using the
%Morlet wavelet. I implemented this using the textbook 'The Illustrated
%Wavelet Transform Handbook' (Paul Addison), pp. 33ff., and the paper
%'Comparison of the Hilbert transform and wavelet methods...' (Le Van
%Quyen, 2001). The parameter 'sigma' is simply the standard deviation of
%the gaussian window (in seconds). The larger this is, the worse the
%temporal resolution, but the better the frequency resolution. If you set sigma=1/(2*pi) seconds,
%then you get a standard deviation in the power spectrum of the Morlet
%wavelet of 1 Hz (but check p. 85 Le Van Quyen paper: they say that
%2/3*(1/sigma) gives frequency resolution).
%'srate' is the sampling rate of the signal (in Hz).'flo' is the lowest frequency
%component we wish to analyze, 'fhi' the highest, and 'deltaf' is the step
%size in frequency (all in Hz).

function [Modulus,Phases,Transform]=morlet_wav(x,srate,sigma,flo,fhi,deltaf)

if(size(x,1)~=1)
    fprintf('Error: Input signal must be a row vector.')
    return
end
N_orig=length(x);
%zero-pad x so that the number of entries is a power of 2, so that the fft
%will be computationally efficient
N=2^(nextpow2(length(x)));
x=[x,zeros(1,N-length(x))];
X=fft(x);

%figure out number of total frequency values at which you will be sampling
%for the time-frequency analysis, and allocate space in 'Transform' (first
%row of 'Transform' contains the power as a function of time for the lowest
%frequency
freqvals=flo:deltaf:fhi;
num_freqvals=length(freqvals);
Transform=zeros(num_freqvals,N);

freq_samples=srate*1/N*(-N/2:N/2-1); %construct array of frequency values at which you sample the Fourier Transform of the wavelet function
for i_f=1:num_freqvals
    
    %construct fourier transform of the Morlet wavelet in such a form that we
    %can use Eq. 2.35 (p. 35, Addison) along with iFFT to determine Transform
    %for specific frequency band. Note that my normalization is not the
    %same as in Addison's textbook.
    W=sqrt(2*pi)*sigma*exp(-2*pi^2*sigma^2*(freq_samples-freqvals(i_f)*ones(1,N)).^2);
    Transform(i_f,:)=ifft(X.*ifftshift(W));
end

%throw away the part of Transform that corresponded to zero-padded portion
%of 'x'
Transform=Transform(:,1:N_orig);

Phases=atan2(imag(Transform),real(Transform));
Modulus=abs(Transform);

return
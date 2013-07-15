function [ output_args ] = heartbeat_preliminary_timing( display, algo, model, time, observ )
%HEARTBEAT_PRELIMINARY_TIMING Get a starting estimate of beat timings

fs = model.fs;

% Filter parameters
fhigh = 2;
flow = 1/3;

% Abs it
observ_abs = abs(observ);

% Filter to a feasible heartbeat range
N = 100;
b = fir1(N, [flow fhigh]/(fs/2));% - ones(1,N+1)/(N+1);
observ_filt = zeros(size(observ_abs));
for ss = 1:4
    observ_filt(ss,:) = filtfilt(b, 1, observ_abs(ss,:)')';
end

% FFT
N = 10000;
f_space = fs/N;
spect = fft(observ_filt', N);
spect(N/2+1:end,:) = [];
freqs = 0:f_space:model.fs/2; freqs(end)=[];

% Find two largest peaks and their phase
mag_spect = sqrt(sum(abs(spect).^2,2));
cond = (mag_spect>[Inf; mag_spect(1:end-1)]) & (mag_spect>[mag_spect(2:end); Inf]);
pi1 = find( mag_spect==max(mag_spect(cond)) ,1);
cond = cond & ((1:N/2)'~=pi1);
pi2 = find( mag_spect==max(mag_spect(cond)) ,1);

freq1 = freqs(pi1);
phase1 = angle(spect(pi1,:));
freq2 = freqs(pi2);
phase2 = angle(spect(pi2,:));

figure, plot(freqs, mag_spect )

end


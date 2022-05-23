%% Test out phase estimation on heart rate signal

load('ECGSignal')
ecg = ECGSignal(:,1);
L = length(ecg);
t = 0:100/L:100-100/L;
fs = 1/mean(diff(t));

L_samp = 5000;
ecg_samp = ecg(1:L_samp);
t_samp = t(1:L_samp);
ecg_filt = lowpass(ecg_samp,10,fs);

figure(1)
clf
plot(t_samp,ecg_filt)
title('ECG Signal')

% Find spectrogram
sig_demean = ecg_filt - mean(ecg_filt);

win_len = round(L_samp/5);
beta = 2;
noverlap = win_len-1;

s = spectrogram(sig_demean,kaiser(win_len,beta),noverlap,1000,fs);
tz = linspace(0,t_samp(L_samp-win_len+1),size(s,2));
fz = linspace(0,fs/2,size(s,1));

figure(2)
clf
mesh(tz,fz(1:100),abs(s(1:100,:)))
view(0,90)


% Obtain phase spectrogram
phase = angle(s);
phase_adj = angle(adj_phase(s,tz,fz));
phase_uw = unwrap2(phase_adj,pi,.2);

figure(3)
clf
mesh(tz,fz(1:100),phase_uw(1:100,:))
view(0,90)

%% Reconstruct signal
testsig = zeros(length(tz),1);
for i = 1:sum(fz<2)
    freqcomp = (abs(s(i,:)).*cos(2*pi*fz(i).*tz+phase_uw(i,:)))';
    testsig = testsig+freqcomp;
end

testsig_filt = lowpass(testsig,1,fs);

ref_sig = sig_demean./max(sig_demean);
recon_sig = testsig./max(testsig_filt);

figure(4)
clf
hold on
plot(t_samp,ref_sig)
plot(tz,recon_sig)

%% Check first frame
frame1 = zeros(win_len,1);
for i = 1:length(fz)
    freqcomp = (abs(s(i,1))*cos(2*pi*fz(i)*t_samp(1:win_len)+phase(i,1)))';
    frame1 = frame1+freqcomp;
end

frame1 = frame1./max(frame1);
ref_sig = sig_demean(1:win_len)./max(sig_demean(1:win_len));

figure(5)
clf
hold on
plot(t_samp(1:win_len),ref_sig)
plot(t_samp(1:win_len),frame1)



%%
% sig1 = sig(1:win_len)-mean(sig(1:win_len));
% sig1_win = sig1 .* kaiser(win_len,beta);
% t1 = t(1:win_len);
% L1 = length(sig1);
% 
% xdft = fft(sig1,L1*2);
% xdft = xdft(1:length(xdft)/2+1);
% xdft = xdft/L1;
% xdft(2:end-1) = 2*xdft(2:end-1);
% fz = linspace(0,fs/2,length(xdft));
% figure(3); plot(fz,abs(xdft));
% xlabel('Hz'); ylabel('Amplitude');
% 
% [s,fz,tz] = czt_spectrogram(
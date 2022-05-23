%% Everything from scratch

% Start with simulated signal (pure cosine)

fs = 1000;
t = 0:1/fs:1-1/fs;
sig = cos(2*pi*10*t);
L = length(sig);

figure(1)
clf
plot(t,sig)

% Find spectrogram
win_len = L/2;
beta = 2;
noverlap = win_len-1;
s = spectrogram(sig,kaiser(win_len,beta),noverlap,fs/2+1);
tz = linspace(0,t_samp(L_samp-win_len+1),size(s,2));
fz = linspace(0,fs/2,size(s,1));

figure(2)
clf
mesh(tz,fz,abs(s))
view(0,90)

% Find unaltered phase
phase = angle(s);

figure(3)
clf
mesh(tz,fz,phase)
view(0,90)

% Adjust phase estimates to compensate for frames
spect_adj = adj_phase2(s,tz,fz);
phase_adj = angle(spect_adj);

figure(4)
clf
mesh(tz,fz,phase_adj)
view(0,90)



function spect_adj = adj_phase2(spect,tz,fz)

dt = mean(diff(tz));
n = 0:length(tz)-1;
spect_adj = spect;

for i = 2:size(spect,1)
    phase_org = angle(spect(i,:));
    amp = abs(spect(i,:));
    nt = 1/(fz(i)*dt)+1;
    phase_shift_uw = 2*pi*n/nt;
    phase_adj = atan(tan(phase_org-phase_shift_uw));
    spect_adj(i,:) = amp.*(cos(phase_adj)+1i*sin(phase_adj));
end
end

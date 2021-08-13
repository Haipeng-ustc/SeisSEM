%% waveform for line AA
clear;clc;close all;

seism_x =  load('seism_data_new_x.dat');
seism_z =  load('seism_data_new_z.dat');
rec_xz  =  load('rec_xz_new.dat');

seism_x = diff(seism_x,1,2);
seism_z = diff(seism_z,1,2);

dt = 2e-3;
nt = 12501;
t  = linspace(0, (nt-1)*dt, nt);
rec_xz  = rec_xz(1,:) - 50000;

seism_x = seism_x(:,1:nt);   seism_x(:,12460:end) = 0.0;
seism_z = seism_z(:,1:nt);   seism_z(:,12460:end) = 0.0;

seism = zeros(2,199,12501);
seism(1,:,:) = seism_x;
seism(2,:,:) = seism_z;


% The station is deployed from 0~247.95km along distance
figure
imagesc(t, rec_xz(1,:), seism_x);  caxis([-1 1]*1e5);
figure
Plot_Trace(rec_xz(:,1:5:end), t(1:nt), seism(:,1:5:end,1:nt), 0.80,'./California_Observation_data_new/waveform.png');

figure 
subplot(211)
plot(t(1:nt), squeeze(seism(1,end,1:nt)), 'k-', 'LineWidth',2.0) ; xlabel('Time (s)'); ylabel('Amplitude');title('velocity-x')
subplot(212)
plot(t(1:nt), squeeze(seism(2,end,1:nt)), 'k-', 'LineWidth',2.0) ; xlabel('Time (s)'); ylabel('Amplitude');title('velocity-z')

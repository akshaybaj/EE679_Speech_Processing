y = wavread('Sounds-u/u_F0-220_F1-300F2-870F3-2240_B1-50B2B3-100_Fs16k.wav');
fs = 16000;
windowing(y,fs,220);
y = wavread('Sounds-u/u_F0-120_F1-300F2-870F3-2240_B1-50B2B3-100_Fs16k.wav');
fs = 16000;
windowing(y,fs,120);
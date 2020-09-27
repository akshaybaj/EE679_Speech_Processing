pkg load signal;

f1=300;
f2=870;
f3=2240;
filename="u";
b0=100;
fs=16000;
f0=[220];

for j=1:columns(f0)
    r1 = exp(-b0*pi*1/fs);
    theta1 = 2*pi*f1*1/fs;
    r2 = exp(-b0*pi*1/fs);
    theta2 = 2*pi*f2*1/fs;
    r3 = exp(-b0*pi*1/fs);
    theta3 = 2*pi*f3*1/fs;
    
    poles1 = [r1*exp(1j*theta1) , r1*exp(-1j*theta1)];
    poles2 = [r2*exp(1j*theta2) , r2*exp(-1j*theta2)];
    poles3 = [r3*exp(1j*theta3) , r3*exp(-1j*theta3)];
    
    b1 = [1 ,0,0];
    a1= [1,-2*r1*cos(theta1),r1**2];
    b2 = [1 ,0,0];
    a2= [1,-2*r2*cos(theta2),r2**2];
    b3 = [1 ,0,0];
    a3= [1,-2*r3*cos(theta3),r3**2];
    b_temp=conv(b1,b2);
    b = conv(b_temp,b3);
    a_temp=conv(a1,a2);
    a = conv(a_temp,a3);
    [h,w] = freqz(b,a) ;
    
    k=figure;
    plot(fs*w/(2*pi),20*log10(abs(h)));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(['Frequency response for filter /',filename,'/']);
    grid on;
    saveas(k,sprintf('output/Frequency_response_%s.png',filename));
    y=input_signal(h,b,a,f0(1,j),fs,0.5,filename);
    windowing(y,w,fs,f0(1,j));
endfor


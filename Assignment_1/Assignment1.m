pkg load signal;

%Question 1 
%f1=900;
%b1=200;
%fs=16000;
%f0=140
%[b,a]=get_coeff(f1,b1,fs);
%h=freq_response(b,a,fs,f1);
%%Question 2
%input_signal(h,b,a,f0,fs,0.5,"900Hz");

%%Question 3(a)
%f0=[120,120,180];
%f1=[300,1200,300];
%b1=[100,200,100];
%for i=1:3
%  [b,a]=get_coeff(f1(1,i),b1(1,i),fs);
%  h=freq_response(b,a,fs,f1(1,i));
%  input_signal(h,b,a,f0(1,i),fs,0.5,num2str(f1(1,i)));
%endfor

%Question 4
f1=[730,270,300];
f2=[1090,2290,870];
f3=[2440,3010,2240];
filename=["a","i","u"];
b0=100;
fs=16000;
f0=[120,220];

for j=1:columns(f0)
  for i=1:columns(f1)
    r1 = exp(-50*pi*1/fs);
    theta1 = 2*pi*f1(1,i)*1/fs;
    r2 = exp(-b0*pi*1/fs);
    theta2 = 2*pi*f2(1,i)*1/fs;
    r3 = exp(-b0*pi*1/fs);
    theta3 = 2*pi*f3(1,i)*1/fs;
    
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
    figure;
    plot(fs*w/(2*pi),20*log10(abs(h)));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(['Frequency response for filter /',filename(1,i),'/']);
    grid on;
    input_signal(h,b,a,f0(1,j),fs,0.5,filename(1,i))
  
  endfor
endfor


function windowing(y,fs,f0)
 m=[0.005,0.010,0.020,0.040];
  
for i=1:columns(m)
  time_=m(1,i);
  M = time_ * fs;
  r=y(600:600+M,1);
  w = .54 + .46*cos(pi*(-M/2:M/2)/M);

  result = r'.*w;
  display(size(result));
  result = fft(result,1024);
  k=figure;
  f = fs*(1:(512))/1024;
  
  plot(f,20*log10(abs(result(1:length(result)/2))));
  xlabel('Frequency(Hz)');
  ylabel('Magnitude (dB)');
  title(['Hamming Window output for F0 =', num2str(f0),' and window length=',num2str(m(1,i)),'s']);
  grid minor;
  saveas(k,sprintf('output/hamming/Hamming_F0_%d_%fs.png',f0,m(1,i)));
  close(k);
endfor

for i=1:columns(m)
  time_=m(1,i);
  M = time_ * fs;
  r=y(600:600+M,1);
  w = ones(1:M,1);
  result = r'.*w;
  display(size(result));
  result = fft(result,1024);
  k=figure;
  f = fs*(1:(512))/1024;
  plot(f,20*log10(abs(result(1:length(result)/2))));
  xlabel('Frequency(Hz)');
  ylabel('Magnitude (dB)');
  title(['Rectangular Window output for F0 =', num2str(f0),' and window length=',num2str(m(1,i)),'s']);
  grid minor;
  saveas(k,sprintf('output/rectangle/Rectangular_F0_%d_%fs.png',f0,m(1,i)));
  close(k);

endfor





function h=freq_response(b,a,fs,f1)
  
  [h,w] = freqz(b,a) ;
  fullname=['assignment1/frequency_response_f1_',num2str(f1),'.jpg']
  figure
  plot(fs*w/(2*pi),20*log10(abs(h)));
  xlabel('Frequency (Hz)');
  ylabel('Magnitude (dB)');
  title(['Frequency response for formant at ',num2str(f1)]);
  grid on;
  
  impulse = zeros(200,1);
  impulse (1,1) = 1;
  y = zeros(200,1);
  [ro,col]=size(y);
  [m,number_of_poles]=size(a);
  for i=number_of_poles:ro
    y(i,1) = b(1,1)*impulse(i-2,1) ;
    for j=2:number_of_poles
      y(i,1)=y(i,1)-a(1,j)*y(i-j+1,1);
    end
  end
  figure;
  fullname=['assignment1/filter_response_f1_',num2str(f1),'.jpg']
  time= linspace(200/fs,1,200);
  plot(time,y);
  xlabel('Time(s)');
  ylabel('Amplitude');
  title(['Impulse Response for formant at ',num2str(f1)] );
  grid on;
  
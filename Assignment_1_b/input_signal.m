function y=input_signal(h,b,a,f0,fs,time_,filename)
  
  t = 0:time_/fs:time_; % 0s to 0.5s with Fs sample freq
  [z,time_length]=size(t);
  x1 = max(0,(square(2*pi*f0*t,0.01)));  % 2*pi* freq * time duration
  %figure;
  %plot(t,x1);
  %title('Triangular pulse train ');
  x1=x1';
  y=zeros(time_length,1);
  
  [ro,col]=size(y);
  [m,number_of_poles]=size(a);
  
  for i=number_of_poles:rows(y)
    y(i,1) = b(1,1)*x1(i-2,1) ;
    for j=2:number_of_poles
      y(i,1)=y(i,1)-a(1,j)*y(i-j+1,1);
    end
  end
  
  
  k=figure;
  t=t';
  plot(t,y);
  xlabel('Time(s)');
  ylabel('Magnitude');
  title(['Filter output for formants corresponding to /',filename,'/ and F0 =', num2str(f0)]);
  grid on;
  saveas(k,sprintf('output/Filter_output_F0_%d.png',f0));
  close(k);
  wavwrite (y,fs, ["output/",filename,"_",num2str(f0),".wav"]);
  

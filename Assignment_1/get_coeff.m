function [b,a]=get_coeff(f1,b1,fs)
  
  r = exp(-b1*pi*1/fs);
  theta = 2*pi*f1*1/fs;
  poles = [r*exp(1j*theta) , r*exp(-1j*theta)];
  b = [1 ,0,0];
  a= [1,-2*r*cos(theta),r**2];
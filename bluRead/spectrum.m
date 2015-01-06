clear all; close all; clc; format compact;

% SET PARAMETERS HERE
  N = 1201;  tmax = 4;
  A0 = 5;  f0 = 9;  phi0 = pi/4;  % tmax*f0 must be an integer for periodicity

% Create t and x vectors
  Nhalf = floor(N/2)+1;
  t = linspace(0,tmax,N);  dt=t(2);
  x = A0*sin(2*pi*f0*t + phi0);

% Fourier transform of x
  xhat = fft(x)/(N);

% two-sided with first entry corresponding to zero frequency
  % -- frequency
  f = linspace(0,1/(2*dt),Nhalf);  df=f(2);  
    if( mod(N,2)==0);  f_all = [f(1:end) -fliplr(f(2:end-1))];
    else;              f_all = [f(1:end) -fliplr(f(2:end))];  end;
  % -- spectrum
  Ex_all = ( xhat.*conj(xhat) ) / df; 
  % -- plot
    figure; plot( f_all, Ex_all ); xlim([-2*f0 2*f0]);
    xlabel('frequency'); ylabel('two-sided power spectral density'); 

% two-sided with center at zero frequency
  % -- frequency
  f_all_c0 = circshift(  f_all, [0 ceil(N/2)-1] );
  % -- spectrum
  Ex_all_c0 = circshift( Ex_all, [0 ceil(N/2)-1] );
  % -- plot
    figure; plot( f_all_c0, Ex_all_c0 ); xlim([-2*f0 2*f0]);
    xlabel('frequency'); ylabel('two-sided power spectral density'); 

% one-sided (double all elements with corresponding negative freq)
  % -- spectrum
  Ex = 2*Ex_all(1:Nhalf);    
    Ex(1) = 0.5*Ex(1);   
    if( mode(N,2)==0); Ex(Nhalf) = 0.5*Ex(Nhalf);  end;
  % -- plot
    figure; plot( f, Ex ); xlim([0 2*f0])
    xlabel('frequency'); ylabel('one-sided power spectral density'); 

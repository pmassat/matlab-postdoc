function [ smoothed, smoothed_D, smoothed_D2 ] = gConvolve( raw, sigma )
%gConvolve smoothes data and derivatives by convoluting with gaussian of SD sigma
%sigma is in number of points

filtergrid=-(10*sigma):(10*sigma);% interval of data smoothing

filter = 1/sqrt(2*pi*sigma^2) * exp(- filtergrid.^2/(2*sigma^2));% gaussian
smoothed=conv(raw,filter,'same');% smoothed data = convolution of data with a gaussian

Dfilter = 1/sqrt(2*pi*sigma^6)*(-filtergrid).* exp(- filtergrid.^2/(2*sigma^2));% first derivative of a gaussian
smoothed_D=conv(raw,Dfilter,'same');% smoothed first derivative of data = convolution of data with derivative of a gaussian

D2filter = 1/sqrt(2*pi*sigma^10)*(filtergrid.^2-sigma^2).* exp(- filtergrid.^2/(2*sigma^2));% second derivative of a gaussian
smoothed_D2=conv(raw,D2filter,'same');% smoothed second derivative of data = convolution of data with second derivative of a gaussian

end


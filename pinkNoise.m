function [ts] = pinkNoise(hurst, ts_length)

n = ts_length;
mu = 0;
sd = 1;
H = hurst;
z = mu + sd.*randn(1, 2*n);
%z = [ 3, 4, 3, 5, 2, 3, 4, 3, 4, 3, 4, 5, 3, 6, 4, 5, 3, 5, 7, 5];
zr = z(1:n);
zi = z(n+1 : 2*n);
zic = -zi;
zi(1) = 0;
zr(1) = zr(1) * sqrt(2);
zi(n) = 0;
zr(n) = zr(n) * sqrt(2);
zr = [zr(1:n), zr(fliplr(2:(n-1)))];
zi = [zi(1:n), zic(fliplr(2:(n-1)))];
z = complex(zr,zi);

k = 0:(n-1);
gammak = ((abs(k-1).^(2*H)) - (2*abs(k).^(2*H)) + (abs(k+1).^(2*H)))/2;
%ind = [0: (n-2), (n-1), fliplr(1:(n-2))];
ind = [0: (n-1), fliplr(1:(n-2))];
gkFGN0 = ifft(gammak(ind + 1) * length(z));
gksqrt = real(gkFGN0);
if( all(gksqrt > 0))
    gksqrt = sqrt(gksqrt);
    z = z.*gksqrt; % matched till here
    z = ifft(z) *length(z);
    z = 0.5*(n-1).^(-0.5)*z;
    ts = real(z);
end

end


% plot(z)
%
% mean_Period =2;
% sd_Period = 0.02;
% for i = 1:length(z)
%     scaled_pink_noise(i) = mean_Period + sd_Period / sd * z(i);
% end
%
% % Can do it without the use of a loop
% scaled_pink_noise1 = mean_Period + (sd_Period/sd) * z;
%
%
% figure(2)
% plot(z)
% hold on;
% plot(scaled_pink_noise);
%
% [x y] = size(scaled_pink_noise)
%
% pink_mean = sum(scaled_pink_noise)/ y
% pink_std = std(scaled_pink_noise)
%
%

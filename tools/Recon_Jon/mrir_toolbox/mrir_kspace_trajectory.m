


N = 16;
n = (1:N)';            % rows

k = sort(rand(1,N)*N); % cols

D = exp(-i*2*pi*(n-1)*(k-1)/N);  % each column vector is a different
                                 % k-space location

figure; imagesc(abs(D'*D/N)); axis image; colorbar;


% TODO: write helper code for reading in gradient waveform dumps from
% IDEA simulator


% TODO: implement conjugate phase reconstruction
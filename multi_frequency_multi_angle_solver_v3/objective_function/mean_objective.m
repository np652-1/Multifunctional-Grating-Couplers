function [fom, gradient] = mean_objective(x,CE,gradients,dirname)
% objective function for multi-angle optimization
% The fom is the sum of coupling efficiencies of all angles
% The gradient is the mean of gradients of all angles 
% CE is a (1, Nth) vector
% gradients is (N, Nth) matrix, N=length(x)
if size(gradients,1)~=length(x)
    error('Error. \ngradient and design variable sizes do not match');
end
fom = mean(CE,[1,2]); % mean over freqs and angles
fprintf('\nfom: %g \n',fom);
if nargout>1
    gradient = mean(gradients,[2,3]); % mean over freqs and angles
end

% Outputs CE, fom, design variable (optional)
% output CE for all pitches
for i=1:size(CE,2)
    fid = fopen(fullfile(dirname,['/CE_pitch_',num2str(i),'.txt']), 'a+');
    fprintf(fid, [repmat('%g\t',1,size(gradients,2)),'\n'], CE(:,i));
    fclose(fid);  
end

% output fom
fid = fopen(fullfile(dirname,'/fom.txt'), 'a+');
fprintf(fid, '%f\n', fom);
fclose(fid);   

% % output structure at cetain condition (# of iterations, or certain fom)
% fid = fopen(fullfile(dirname,'/fom.txt'), 'r');
% foms = fscanf(fid, '%f');
% fclose(fid);
% iter = length(foms); % number of iterations
% if mod(iter,10)==0
%     % output optimized structure
%     writematrix(x(:), fullfile(dirname,['gc_density_iter_',num2str(iter,'%d'),'.txt']));
% end
end
function cmap = singleColorIntensityColormap(wavelength)
    % Creates an intensity colormap for a given wavelength
    % Zero intensity is shown in white, Stronger intensity is shown in
    % stronger colors
    % Input: wavelength, unit: nm
    % Output: colormap used by MATLAB
    
    % Define the number of levels in the colormap
    nLevels = 256;

    % Define the target color
    targetColor = wavelengthToRGB(wavelength);

    % Generate the colormap
    cmap = zeros(nLevels, 3); % Preallocate
    for i = 1:nLevels
        intensity = (i-1) / (nLevels-1); % Scale intensity from 0 to 1
        cmap(i, :) = (1-intensity) * [1, 1, 1] + intensity * targetColor; % Interpolate between white and the target color
    end
end


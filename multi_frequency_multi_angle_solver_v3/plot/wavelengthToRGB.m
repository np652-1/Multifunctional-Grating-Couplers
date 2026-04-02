function rgb = wavelengthToRGB(wavelength)
    % wavelengthToRGB converts a wavelength in nanometers to an RGB triplet.
    %
    % Input:
    %   wavelength - Wavelength in nanometers (should be between 380 and 780)
    %
    % Output:
    %   rgb - 1x3 vector representing the RGB color (values in the range [0, 1])
    % This function is written by ChatGPT with hint:
    % Write a function that returns rbg values of a wavelength
    % This may not be accurate and can be tuned if necessary.

    if wavelength < 380 || wavelength > 780
        error('Wavelength must be in the range 380 to 780 nm.');
    end

    % Initialize RGB components
    R = 0;
    G = 0;
    B = 0;

    % Map wavelength to RGB based on visible spectrum
    if wavelength >= 380 && wavelength < 440
        R = -(wavelength - 440) / (440 - 380);
        B = 1;
    elseif wavelength >= 440 && wavelength < 490
        G = (wavelength - 440) / (490 - 440);
        B = 1;
    elseif wavelength >= 490 && wavelength < 510
        G = 1;
        B = -(wavelength - 510) / (510 - 490);
    elseif wavelength >= 510 && wavelength < 580
        R = (wavelength - 510) / (580 - 510);
        G = 1;
    elseif wavelength >= 580 && wavelength < 645
        R = 1;
        G = -(wavelength - 645) / (645 - 580);
    elseif wavelength >= 645 && wavelength <= 780
        R = 1;
    end

    % Apply intensity correction for low/high wavelengths
    if wavelength >= 380 && wavelength < 420
        alpha = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);
    elseif wavelength > 700 && wavelength <= 780
        alpha = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);
    else
        alpha = 1;
    end

    % Combine RGB components and apply intensity correction
    rgb = [R, G, B] * alpha;

    % Ensure RGB values are within [0, 1]
    rgb = max(0, min(1, rgb));
end

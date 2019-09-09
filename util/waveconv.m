function [wave_cos, wave_sin] = waveconv(bin_resolution, virtual_wavelength, cycles, data)
    c = 3e8;
    s_z = bin_resolution * c;
    samples = round(cycles * virtual_wavelength / (bin_resolution * c));
    num_cycles = samples * s_z / virtual_wavelength;
    sigma = 0.3;
    
    % generate sin/cos signals
    sin_wave = sin(2*pi*(num_cycles * linspace(1,samples,samples)')/samples);
    cos_wave = cos(2*pi*(num_cycles * linspace(1,samples,samples)')/samples);
    
    % apply window function
    window = single(gausswin(samples, 1/sigma));  
    virtual_wave_sin = sin_wave .* window;
    virtual_wave_cos = cos_wave .* window;

    % convolve with measurements
    wave_cos = single(zeros(size(data)));
    wave_sin = single(zeros(size(data)));
    for laser_index = 1 : size(data,2)
        for camera_index = 1 : size(data,3)
            time_response = squeeze(data(:,laser_index, camera_index));
            tmp_real = conv(time_response,virtual_wave_sin, 'same');
            tmp_img = conv(time_response,virtual_wave_cos, 'same');
            wave_cos(:,laser_index, camera_index) = tmp_real;
            wave_sin(:,laser_index, camera_index) = tmp_img;
        end
    end 
end






function [real, img] = waveconvolution(time_resolution, virtual_guassian_size, virtual_lamda, data_t)
% Input parameter: 
% time_resolution:  unit second, time resolution per time bin 
% virtual_guassian_size:   integer, represent as virtual phasor wave guassian envelop
% virtual_lamda:     unit cm, virtual phasor wavelength
% data_t:               3D temporal measurement

% Output parameter:
% real, img: real and imaginary part of the virtual phasor wave

    c = 3e8; % speed of light
    s_z = time_resolution * c * 100; % unit: cm

    Sinusoid_pattern = virtual_guassian_size * s_z / virtual_lamda;

    Gauss_sigma = 0.3;

    Sinusoid_pattern = single(Sinusoid_pattern);
    virtual_guassian_size = single(virtual_guassian_size);
    
    sin_wave = sin(2*pi*(Sinusoid_pattern * linspace(1,virtual_guassian_size,virtual_guassian_size)')/virtual_guassian_size);
    cos_wave = cos(2*pi*(Sinusoid_pattern * linspace(1,virtual_guassian_size,virtual_guassian_size)')/virtual_guassian_size);
    gauss_wave = single(gausswin(virtual_guassian_size, 1/Gauss_sigma));

    Virtual_Wave_sin = sin_wave .* gauss_wave;
    Virtual_Wave_cos = cos_wave .* gauss_wave;

    real = single(zeros(size(data_t)));
    img = single(zeros(size(data_t)));
    for laser_index = 1 : size(data_t,2)
        for camera_index = 1 : size(data_t,3)
            time_response = squeeze(data_t(:,laser_index, camera_index));

            % Wave convolution
            tmp_real = conv(time_response,Virtual_Wave_sin, 'same');
            tmp_img = conv(time_response,Virtual_Wave_cos, 'same');
            real(:,laser_index, camera_index) = tmp_real;
            img(:,laser_index, camera_index) = tmp_img;

        end
    end 

end

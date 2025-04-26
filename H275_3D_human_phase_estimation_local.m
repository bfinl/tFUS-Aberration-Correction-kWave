%%%
% SPDX-License-Identifier: LGPL-3.0-or-later
% Copyright (c) 2025 Zherui Li
%%%
%%

clearvars;

addpath(genpath('./k-wave-toolbox-version-1.4'))
addpath(genpath('./H-275 3D Simulations'))

set(0,'DefaultFigureVisible', 'off')
tStart = tic;

%%
% =========================================================================
% DEFINE LITERALS
% =========================================================================
    
% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: MATLAB GPU code
%   3: C++ code
%   4: CUDA code
%   5: Save H5 file only
model           = 4;

% subject ID
subjectID = 'tFUS47';

% medium parameters
c0              = 1482;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]
alpha0          = 1/11.5128*3.48e-4;  % power law absorption coefficient [dB/(MHz^y cm)]

% skull parameters
c_skull         = 2850;
rho_skull       = 1732;
alpha_skull     = 1/11.5128*85;

% source parameters
source_f0       = 7e5;      % source frequency [Hz]
source_diameter = 4e-3;     % disc diameter [m]
source_amp      = 1e6;      % source pressure [Pa]
focus_depth     = 35e-3;    % focus depth (ROC of the transducer) [m]

% grid parameters
axial_size      = 70e-3;    % total grid size in the axial dimension [m]
lateral_size    = 150e-3;   % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 2;        % number of points per wavelength
t_end           = 50e-6;    % total compute time [s]
record_periods  = 1;        % number of periods to record
cfl             = 0.02;     % CFL number
source_x_offset = 20;       % grid points to offset the source
bli_tolerance   = 0.01;     % tolerance for truncation of the off-grid source points
upsampling_rate = 10;       % density of integration points relative to grid
%%
% =========================================================================
% DEFINE KGRID
% =========================================================================

% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * source_f0);   % [m]

% compute the size of the grid
Nx = roundEven(axial_size / dx) + source_x_offset;
Ny = roundEven(lateral_size / dx);
Nz = Ny;

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% compute points per temporal period
PPP = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (PPP * source_f0);

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);
%%
% =========================================================================
% DEFINE KWAVEARRAY FOR TRANDUCER
% =========================================================================

% Load H-275 element sptial positions
pos = load('./H-275 3D Simulations/H275_elementsPos_Sorted.mat');
pos = getfield(pos, 'H275_ElementPos_sorted'); % #, X, Y, Z
pos(:,[2 4]) = pos(:,[4 2]); % #, Z, Y, X
pos(:,[3 4]) = pos(:,[4 3]); % #, Z, X, Y

% set element position
elem_pos = pos(:,2:4) * 1e-3;
elem_pos(:,1) = elem_pos(:,1) + abs(min(gather(elem_pos(:,1)))) + ... % let all x>=0
    kgrid.x_vec(1) + (source_x_offset-15) * kgrid.dx;

% set element focus position
focus_pos = [min(gather(elem_pos(:,1))) + focus_depth, 0, 0];

% create empty kWaveArray
karray = kWaveArray('BLITolerance', bli_tolerance, 'UpsamplingRate', upsampling_rate, 'SinglePrecision', true);

% add disc shaped elements
for i = 1:size(elem_pos,1)
    karray.addDiscElement(elem_pos(i,:), source_diameter, focus_pos);
end
%%
% =========================================================================
% DEFINE BINARY SKULL MASK
% =========================================================================

% load skull model from pCT
skull_bone = niftiread(['./H-275 3D Simulations/subjects/' subjectID '/' subjectID '_pct.nii']);

%%% cut skull bone
% V5
skull_bone(88:end, :, :) = 0;
skull_bone(:, 100:end, :) = 0;
skull_bone(:, :, 1:70) = 0;

skull_bone = skull_bone > 140; % set threshold for binary mask

scale_factor = 1 / (dx*1e3);
skull_bone = imresize3(skull_bone, scale_factor, 'nearest');

%%% set valid range
% V5
skull_bone = skull_bone(7:81, :, :); % x-axis
skull_bone = skull_bone(:, 24:92, :);  % y-axis
skull_bone = skull_bone(:, :, 65:215); % z-axis

%%% rotation
% V5
skull_bone = imrotate3(double(skull_bone), -90, [0 0 1], 'nearest');
skull_bone = imrotate3(double(skull_bone), 45, [0 0 1], 'nearest');

skull_mask = zeros(Nx, Ny, Nz, 'logical');
skull_mask(ceil(Nx/2)-ceil(size(skull_bone,1)/2)+34:ceil(Nx/2)+ceil(size(skull_bone,1)/2)-46+34, ...
    ceil(Ny/2)-ceil(size(skull_bone,2)/2)+0:ceil(Ny/2)+ceil(size(skull_bone,2)/2)-1+0, ...
    ceil(Nz/2)-ceil(size(skull_bone,3)/2)+16:ceil(Nz/2)+ceil(size(skull_bone,3)/2)-24+16) = skull_bone(29:85,:,1:end-22);
%%
% =========================================================================
% SIMULATIONS
% =========================================================================

% create matrix to save delay number
data_to_write = zeros(size(pos,1), 3);
data_to_write(:,1) = pos(:,1); % 1st value in each row: element index

% write NaN to the existing phase file only if it does not exist
filename = ['./H-275 3D Simulations/subjects/' subjectID '/' subjectID '_phase_delay_v5.xlsx'];

if ~exist(filename, 'file')
    hearder_to_write = ["Element#","TimeStepDiff","PhaseDelay"];
    nan_to_write = repmat("NaN",128,3);

    writematrix(hearder_to_write, filename, 'Sheet', 1, 'Range', 'A1:C1');
    writematrix(nan_to_write, filename, 'Sheet', 1, 'Range', 'A2:C129');
    fprintf('Phase file created and NaN values written: %s\n', filename);
    disp('Phase estimation loop begins ...');
else
    fprintf('Phase file already exists: %s\nContinue the simulation that may have been interrupted ...\n', filename);
end

% define path to store figures and measured single element sensor data
figureFolderPath = ['./H-275 3D Simulations/figs/human/subjects/' subjectID '/'];
sensorFolderPath = ['./H-275 3D Simulations/saved workspace/human/subjects/' subjectID '/'];
% define path to load source struct
sourceFolderPath = './H-275 3D Simulations/subjects/source/';

% check if the folders exist
if ~exist(figureFolderPath, 'dir')
    % if it doesn't exist, create it
    mkdir(figureFolderPath);
    fprintf('Figure folder created: %s\n', figureFolderPath);
else
    fprintf('Figure folder already exists: %s\n', figureFolderPath);
end

if ~exist(sensorFolderPath, 'dir')
    mkdir(sensorFolderPath);
    fprintf('Sensor data folder created: %s\n', sensorFolderPath);
else
    fprintf('Sensor data folder already exists: %s\n', sensorFolderPath);
end

%%
% =========================================================================
% LOOP TO OBATIN ALL PHASE DELAY NUMBERS
% =========================================================================

for i = 1:size(pos,1) % iterate all elements
    % get current element index in case they are not sorted
    elem_ind = pos(i,1);

    % check if phase estimation has once been done for current element
    phase_data_in_file = readmatrix(filename,'Sheet',1,'Range','A2:C129');
    if ~isnan(phase_data_in_file(i,1))
        disp(['Simulation for Element #' num2str(elem_ind) ' has already been done (' num2str(i) '/' num2str(size(pos,1)) ')!']);
        continue;
    end

    % simulate w/ and w/o skull to obtain 2 waveforms at target/focus
    for j = 1:2 % define if w/(2) or w/o(1) skull
        %%
        % =========================================================================
        % DEFINE SOURCE SIGNAL
        % =========================================================================
        
        % if j==1 % create source only once
        %     % create time varying source
        %     source_sig = createCWSignals(kgrid.t_array, source_f0, source_amp, 0);
        %     source_sig_all = repmat(source_sig, size(elem_pos,1), 1); % all elements
        %     % only activate single element once at a time
        %     if i==1
        %         source_sig_all(2:end, :) = 0;
        %     elseif i==size(pos,1)
        %         source_sig_all(1:end-1, :) = 0;
        %     else
        %         source_sig_all(1:i-1, :) = 0;
        %         source_sig_all(i+1:end, :) = 0;
        %     end
        %     
        %     % assign binary mask
        %     source.p_mask = karray.getArrayBinaryMask(kgrid);
        %     
        %     % assign source signals
        %     source.p = karray.getDistributedSourceSignal(kgrid, source_sig_all);
        % end

        if j==1 % load source only once
            load(string([sourceFolderPath 'H275_source_single_' num2str(elem_ind) '.mat']));
            disp(['Source loaded for Element #' num2str(elem_ind)]);
        end
        %%
        % =========================================================================
        % DEFINE MEDIUM
        % =========================================================================
        
        % assign medium properties
        if j==1 % without skull
            medium.sound_speed = c0.*ones(Nx, Ny, Nz);
            medium.density = rho0.*ones(Nx, Ny, Nz);
            medium.alpha_coeff = alpha0.*ones(Nx, Ny, Nz);
        else % with skull
            medium.sound_speed = c0.*ones(Nx, Ny, Nz);
            medium.sound_speed(skull_mask) = c_skull;
            medium.density = rho0.*ones(Nx, Ny, Nz);
            medium.density(skull_mask) = rho_skull;
            medium.alpha_coeff = alpha0.*ones(Nx, Ny, Nz);
            medium.alpha_coeff(skull_mask) = alpha_skull;
        end
        
        medium.alpha_power = 1.001;
        %%
        % =========================================================================
        % DEFINE SENSOR
        % =========================================================================
        
        % set sensor mask to record central plane
        sensor.mask = zeros(Nx, Ny, Nz);
        sensor.mask(source_x_offset + 2:end, :, Nz/2 + 1) = 1;
        
        % record the pressure
        sensor.record = {'p'};
        
        sensor.record_start_index = 1; % record all time pts
        %%
        % =========================================================================
        % SIMULATION
        % =========================================================================
        
        % set input options
        input_args = {...
            'PMLSize', 'auto', ...
            'PMLInside', false, ...
            'PlotPML', false, ...
            'DisplayMask', 'off', ...
            'PlotSim', false};
        
        % run code
        switch model
            case 1
                
                % MATLAB CPU
                sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                    input_args{:}, ...
                    'DataCast', 'single', ...
                    'PlotScale', [-1, 1] * source_amp);
                
            case 2
                
                % MATLAB GPU
                sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                    input_args{:}, ...
                    'DataCast', 'gpuArray-single', ...
                    'PlotScale', [-1, 1] * source_amp);
                
            case 3
                
                % C++
                sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
                
            case 4
                
                % C++/CUDA GPU
                sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
                
            case 5
        
                % Save H5 file only
                filename = './test_input.h5';
                sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                    'SaveToDisk', filename);
        end
        %%
        % extract amplitude from the sensor data
        amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
            'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
        
        % reshape data
        amp = reshape(amp, Nx - source_x_offset - 1, Ny);
        
        % define axis vectors for plotting
        x_vec = kgrid.x_vec(source_x_offset + 2:end, :) - kgrid.x_vec(source_x_offset + 1);
        y_vec = kgrid.y_vec;
        %%
        % =========================================================================
        % VISUALISATION
        % =========================================================================
        
        % plot the pressure field
        if j==1 % without skull
            figure;
            imagesc(1e3 * y_vec, 1e3 * x_vec, amp);
            c = colorbar;
            c.Label.String = "Pressure (Pa)";
            xlabel('Lateral Position [z, mm]');
            ylabel('Axial Position [x, mm]');
            axis image;
            title(['Free-Field Pressure Field (Element #' num2str(elem_ind) ')']);
            
            f = gcf;
            exportgraphics(f,[figureFolderPath 'P_free-field_single_' num2str(elem_ind) '.png'],'Resolution',300);
        else % with skull
            skull = imgaussfilt(squeeze(medium.density(source_x_offset + 2:end, :, Nz/2 + 1)), 0.8);
            
            figure
            overlayPlot(1e3 * y_vec, 1e3 * x_vec, skull, amp, 'LogComp', false, ...
                'PlotScale', [0 gather(max(amp(:)))], ...
                'ColorBar', true, 'ColorBarTitle', "Pressure (Pa)");
            axis image;
            colormap(jet);
            xlabel('Lateral Position [z, mm]');
            ylabel('Axial Position [x, mm]');
            title(['Overlayed Pressure Field (Element #' num2str(elem_ind) ')']);

            f = gcf;
            exportgraphics(f,[figureFolderPath 'P_skull_single_' num2str(elem_ind) '.png'],'Resolution',300);
        end
        %%
        % =========================================================================
        % WAVEFORM MEASUREMENT
        % =========================================================================

        % get pressure at focus vs time
        p_focus = zeros(1,size(sensor_data.p,2));
        for k = 1:size(sensor_data.p,2)
            p_tmp = sensor_data.p(:,k);
            p_tmp = reshape(p_tmp, Nx - source_x_offset - 1, Ny);
            p_focus(k) = p_tmp(17,71);
        end
        
        if j==1 % without skull
            p_free_field = p_focus;
            save(string([sensorFolderPath 'H275_P_T_free-field_single_' num2str(elem_ind) '.mat']), "p_free_field", "-v7.3");
        else % with skull
            p_skull = p_focus;
            save(string([sensorFolderPath 'H275_P_T_skull_single_' num2str(elem_ind) '.mat']), "p_skull", "-v7.3");
        end
    end
    
    %%
    % =========================================================================
    % PHASE DELAY ESTIMATION
    % =========================================================================
    
    % nomarlize pressure data to [-1 1]
    p_min_single = min(p_free_field);
    p_max_single = max(p_free_field);
    p_max_skull_single = max(p_skull);
    norm_p_free_field = 2*(p_free_field-p_min_single)/(p_max_single-p_min_single)-1;
    norm_p_skull = 2*(p_skull-p_min_single)/(p_max_single-p_min_single)-1;

    % find peaks
    [maxv_free_field, maxl_free_field] = findpeaks(norm_p_free_field, ...
        'minpeakdistance',10, ...
        'MinPeakHeight',0.005);
    
    [maxv_skull, maxl_skull] = findpeaks(norm_p_skull, ...
        'minpeakdistance',10, ...
        'MinPeakHeight',0.005);
    
    % calculate phase delay
    if size(maxl_skull,2)>0 && size(maxl_free_field,2)>0
        time_step_diff = maxl_skull(end) - maxl_free_field(end);
        phase_delay = 2 * pi * source_f0 * time_step_diff * kgrid.dt; % phi=2pi*freq*time_delay
        % plot waveform compare
        figure
        plot(norm_p_free_field);
        hold on
        plot(norm_p_skull);
        hold on
        plot(maxl_free_field(end), maxv_free_field(end), '*', 'color', 'R');
        hold on
        plot(maxl_skull(end), maxv_skull(end), '*', 'color', 'G');
        hold off
        xlim([0 size(norm_p_free_field,2)]);
        xlabel('Time Step');
        ylabel('Normalized Pressure');
        title(['Pressure At Focus (Element #' num2str(elem_ind) ')']);
        legend(["Free-Field", "With Skull"],'Location','northwest');
    else % cannot find peaks in the condition with skull
        % set both time diff and phase delay to -1 as a flag
        time_step_diff = -1;
        phase_delay = -1;
        % plot waveform compare
        figure
        plot(norm_p_free_field);
        hold on
        plot(norm_p_skull);
        if size(maxl_free_field,2)>0
            hold on
            plot(maxl_free_field(end), maxv_free_field(end), '*', 'color', 'R');
        end
        hold off
        xlim([0 size(norm_p_free_field,2)]);
        xlabel('Time Step');
        ylabel('Normalized Pressure');
        title(['Pressure At Focus (Element #' num2str(elem_ind) ')']);
        legend(["Free-Field", "With Skull"],'Location','northwest');
    end

    data_to_write(i,2) = time_step_diff; % 2nd value in each row: time diff
    data_to_write(i,3) = phase_delay; % 3rd value in each row: phase delay
    
    f = gcf;
    exportgraphics(f,[figureFolderPath 'P_compare_single_' num2str(elem_ind) '.png'],'Resolution',300);

    % write data for current element to file
    writematrix(data_to_write(i,:),filename,'Sheet',1,'Range',['A' num2str(i+1) ':C' num2str(i+1)]);
    disp(['Phase data saved for Element #' num2str(elem_ind)]);

    disp(['Simulation for Element #' num2str(elem_ind) ' has finished (' num2str(i) '/' num2str(size(pos,1)) ') ...']);
end

disp(['Simulation All Finished! Phase data all saved to: ' filename]);

% display total simulation time
tEnd = toc(tStart);
disp(['Total Simulation Time: ',second_change(tEnd)]);

% set figure back to visible
set(0,'DefaultFigureVisible', 'on')

%%
% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

% Convert seconds to HH:MM:SS string
function Output = second_change(num)
    hour = floor(num/3600);
    minute = floor(mod(num,3600)/60);
    second = num - 3600*hour - 60*minute;
    
    if hour < 10
        hour = ['0',mat2str(hour)];
    else
        hour = mat2str(hour);
    end
    
    if minute < 10
        minute = ['0',mat2str(minute)];
    else
        minute = mat2str(minute);
    end
    
    if second < 10
        second = roundn(second,-2);
        second = ['0',mat2str(second)];
    else
        second = roundn(second,-2);
        second = mat2str(second);
    end
    
    Output = [hour,':',minute,':',second];
end
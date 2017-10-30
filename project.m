%% Settings

steps = 700;
steps_per_sec = 60; % steps per second in video
fps = 30; % frames per second in video
max_fertile_dist = 51; % max distance from sea before 0 fertility (change this for different images)
specify_points = true; % whether to specify initial clusters manually
load_fertility = true; % whether to load fertility from file
map_noise = 10; % random noise in base map image
num_clusters = 40; % number of cells in a starting cluster

% behaviour conditions
ideal_neighbours = 0; % survival chance goes down with every deviation from this number
% fertility at which survival will be multiplied by 1
% values above this will increase survival chance
fertility_point = 0.75;
temperature_point = 0.45; % temperature at and above which survival will be multiplied by 1
movement_chance = 0.4; % chance per step that a cell will move
landing_chance = 0.5; % chance of ship landing when encountering land
voyage_base_chance = 0.001; % initial chance of setting on a voyage when near water
voyage_max_chance = 0.015; % maximium chance of above
change_course_chance = 0.1; % chance of changing course mid voyage

% display/output
green = [107 238 67]; % colour of fertile land
brown = [131 80 15]; % colour of barren land
sea = [25 95 220]; % colour of sea
heat_value = 2; % red added per step to heat map when a cell is over a position in the grid 
bw_output = true; % outputs black/white map with heat
colour_output = false; % outputs colour map
heat_output = false; % outputs black and white with heat without cells
pop_output = false; % outputs population graph over time
sound_output = true; % outputs sine wave sound based on population growth
sound_fs = 8192;

%% Generate Grid

% load image
orig_image = imread('earth.png');
grid = double(orig_image);
grid_size = size(grid);

% this just uses red channel for now; need grayscale with typecasting out of
% uint8 to work for other images
grid(:, :, 1) = grid(:, :, 1) < 120; % land/sea
grid(:, :, 2) = 0; % cell present - 0 for now

half_height = grid_size(1) / 2;
sea_borders = [];

if load_fertility
	load('fertility.mat');
	grid(:, :, 3) = fertility;
else
	% find every sea pixel with bordering land
	for y = 1:grid_size(1)
    	for x = 1:grid_size(2)
        	if ~grid(y, x, 1)
            	tiles = grid(max(y - 1, 1):min(y + 1, grid_size(1)), max(x - 1, 1):min(x + 1, grid_size(2)), 1);
            	if any(tiles)
                	sea_borders(end + 1, :) = [x y];
            	end
        	end
    	end
	end
    
	% this code takes a *long* time to run
	% add option later on to use precomputed data
	for y = 1:grid_size(1)
    	for x = 1:grid_size(2)
        	if grid(y, x, 1)
            	% compute distance with sea border tile and take the minimum
            	dist = sqrt((sea_borders(:, 1) - x).^2 + (sea_borders(:, 2) - y).^2);
            	grid(y, x, 3) = max(min(1 - min(dist) / max_fertile_dist, 1), 0);
        	end
    	end
	end
end


y_mat = repmat((1:grid_size(1))', 1, grid_size(2));
x_mat = repmat(1:grid_size(2), grid_size(1), 1);

% temperature relative to equator (via y axis)
grid(:, :, 4) = max(min(1 - abs(half_height - y_mat)*1.15 / half_height, 1), 0);

% survival factor (combines fertitlity and temperature)
grid(:, :, 5) = (0.984 + 0.016 * min(grid(:, :, 3) / fertility_point, 1)) ...
              .* (0.9775 + 0.0215 * min(grid(:, :, 4) / temperature_point, 1)) ...
              .* grid(:, :, 1);

grid(:, :, 6) = x_mat; % reference to x position
grid(:, :, 7) = y_mat; % reference to y position

%% Generate Base Image

base_image = uint8(zeros(grid_size));
fertility_factor = max(grid(:, :, 3) - 0.15, 0);

% land/sea, fertility, noise
% formula: is_land * fertility_gradient_between_brown_and_green ...
%      	+ is_sea * sea_colour
base_image(:, :, 1) = grid(:, :, 1) .* (brown(1) + (green(1) - brown(1)) * fertility_factor) ...
                    + ~grid(:, :, 1) * sea(1);
base_image(:, :, 2) = grid(:, :, 1) .* (brown(2) + (green(2) - brown(2)) * fertility_factor) ...
                    + ~grid(:, :, 1) * sea(2);
base_image(:, :, 3) = grid(:, :, 1) .* (brown(3) + (green(3) - brown(3)) * fertility_factor) ...
                    + ~grid(:, :, 1) * sea(3);

% temperature
% combines original land with grayscale sum according to a temperature factor
temp_factor = max((0.6 - grid(:, :, 4)) / 0.6, 0) .* grid(:, :, 1);
gray_sum = double(base_image(:, :, 1)) + double(base_image(:, :, 2)) + double(base_image(:, :, 3));
base_image(:, :, 1) = double(base_image(:, :, 1)) .* (1 - temp_factor) + gray_sum .* temp_factor;
base_image(:, :, 2) = double(base_image(:, :, 2)) .* (1 - temp_factor) + gray_sum .* temp_factor;
base_image(:, :, 3) = double(base_image(:, :, 3)) .* (1 - temp_factor) + gray_sum .* temp_factor;

% random noise
base_image = base_image + uint8(-map_noise / 2 + rand(size(base_image)) * map_noise);

%% Populate Initial Clusters

% cells should be stored in a linked list (working on making one)
cells = CellList;

if specify_points
    figure
    imshow(base_image)
    [cx, cy] = getpts;
    
    for i = 1:length(cx)
        for y = (cy(i) - 1):(cy(i) + 1)
            for x = (cx(i) - 1):(cx(i) + 1)
                y = floor(y);
                x = floor(x);
                
                if grid(y, x, 1)
                   cells.add(x, y);
                   grid(y, x, 2) = 1;
                end
            end
        end
    end
else
    % cluster initial positions
    for i = 1:num_clusters
        guess_pos = [ceil(rand*grid_size(2)) ceil(rand*grid_size(1))];

        while guess_pos(1) == grid_size(2) || guess_pos(1) == 1 || guess_pos(2) == grid_size(1) ...
            || guess_pos(2) == 1 || grid(guess_pos(2), guess_pos(1), 5) < 0.5
            guess_pos = [ceil(rand*grid_size(2)) ceil(rand*grid_size(1))];
        end

        for y = (guess_pos(2) - 1):(guess_pos(2) + 1)
            for x = (guess_pos(1) - 1):(guess_pos(1) + 1)
                if grid(y, x, 1)
                   cells.add(x, y);
                   grid(y, x, 2) = 1;
                end
            end
        end
    end
end

%% Create Video

% initialise video
if bw_output
    video1 = VideoWriter('sim_bw');
    video1.FrameRate = fps;
    open(video1);
end

if colour_output
    video2 = VideoWriter('sim_colour');
    video2.FrameRate = fps;
    open(video2);
end

if heat_output
    video3 = VideoWriter('sim_heat');
    video3.FrameRate = fps;
    open(video3);
end

if pop_output
    video4 = VideoWriter('sim_population');
    video4.FrameRate = fps;
    open(video4);
end

render_interval = round(steps_per_sec / fps);

if sound_output
    bits_per_step = round(sound_fs / steps_per_sec);
    sound = zeros(bits_per_step * steps, 1);
end

% neighbouring indicies
% wraps around
north = [grid_size(1) 1:(grid_size(1) - 1)];
south = [2:grid_size(1) 1];
east  = [2:grid_size(2) 1];
west  = [grid_size(2) 1:(grid_size(2) - 1)];

% cell counts at each step
population = zeros(1, steps);

for i = 1:steps
	% sample loop code using the linked list
	c = cells.First;
    population(i) = cells.Length;
    add_cells = [];
    
	while ~islogical(c)
        % with hard boundaries:
        % grid(north(c.Y):south(c.Y), west(c.X):east(c.X), :);
        land = [
            grid(north(c.Y), west(c.X), :) grid(north(c.Y), c.X, :) grid(north(c.Y), east(c.X), :);
            grid(c.Y, west(c.X), :) grid(c.Y, c.X, :) grid(c.Y, east(c.X), :);
            grid(south(c.Y), west(c.X), :) grid(south(c.Y), c.X, :) grid(south(c.Y), east(c.X), :)
        ];
        
        neighbour_count = sum(sum(land(:, :, 2))) - 1;
        
        if land(2, 2, 1) == 1
            survival_chance = (1 - abs(neighbour_count - ideal_neighbours) * 0.003) * (0.5 + 0.5 * land(2, 2, 5));
        else
            survival_chance = 1 - neighbour_count * 0.001 - (1 - min(land(2, 2, 4) / temperature_point, 1)) * 0.02;
        end
        
        if rand > survival_chance
            cells.remove(c);
            grid(c.Y, c.X, 2) = 0;
            c = c.Next;
            continue
        end
        
        empty_land = land(:, :, 1) .* ~land(:, :, 2);
        empty_land(:, :, 2:4) = land(:, :, 5:7);
        empty_sea = ~land(:, :, 1) .* ~land(:, :, 2);
        empty_sea(:, :, 2:3) = land(:, :, 6:7);
        
        if c.PX == c.X && c.PY == c.Y
            c.Times = c.Times + 1;
        end
        
        if land(2, 2, 1) == 1
            % no reproduction or movement if all immediate spaces are full
            if neighbour_count < 8
                % chance of voyage increases with time
                if c.CanLand == 0 && rand < min(voyage_base_chance * (i / 350), voyage_max_chance) ...
                   && any(any(empty_sea(:, :, 1)))
                    [y, x] = find(empty_sea(:, :, 1));
                    index = floor(rand * (length(y) - 1) + 1);
                    yi = y(index);
                    xi = x(index);
                    grid(c.Y, c.X, 2) = 0;
                    c.X = empty_sea(yi, xi, 2);
                    c.Y = empty_sea(yi, xi, 3);
                    grid(c.Y, c.X, 2) = 1;
                    c.DepartX = c.X;
                    c.DepartY = c.Y;
                    c.Course = atan2(-(yi - 2), xi - 2) - pi / 6 + pi / 3 * rand;
                    c.CanLand = 4;
                else
                    % determine empty land cells
                    empty_cells = [];

                    for y = 1:3
                        for x = 1:3
                            if empty_land(y, x, 1)
                                empty_cells(end + 1, :) = empty_land(y, x, 2:4);
                            end
                        end
                    end

                    num_empty_cells = size(empty_cells);
                    num_empty_cells = num_empty_cells(1);

                    % reproduction
                    if num_empty_cells
                        reproduction_chance = .0052 * (1 + neighbour_count * 0.013);      

                        if rand < reproduction_chance
                           spawn_i = floor(rand * (num_empty_cells - 1) + 1);
                           pos = empty_cells(spawn_i, 2:3);
                           grid(pos(2), pos(1), 2) = 1;
                           empty_cells(spawn_i, :) = [];
                           num_empty_cells = num_empty_cells - 1;
                           add_cells(end + 1, :) = pos;
                        end
                    end

                    % movement
                    if num_empty_cells % check if reproduction changed anything
                        if rand < movement_chance
                            min_survival = min(empty_cells(:, 1));
                            survival_range = max(empty_cells(:, 1)) - min_survival;

                            r = rand;
                            base_chance = 1 / num_empty_cells;
                            chance = 0;

                            for j = 1:num_empty_cells
                                % part normal chance distribution, part weighted by
                                % range of survival factors in neighbours
                                chance = chance + 0.7 * base_chance + 0.3 * base_chance ...
                                       * (1 - (survival_range - (empty_cells(j, 1) - min_survival)));

                                if r < chance
                                    pos = empty_cells(j, 2:3);
                                    grid(pos(2), pos(1), 2) = 1;
                                    grid(c.Y, c.X, 2) = 0;
                                    c.X = pos(1);
                                    c.Y = pos(2);
                                    break
                                end
                            end     
                        end
                    end
                end
            end
        else
            landed = false;
            
            % chance of landing
            if c.CanLand == 0 && any(any(empty_land(:, :, 1)))
                if rand < landing_chance
                    [y, x] = find(empty_land(:, :, 1));
                    index = floor(rand * (length(y) - 1) + 1);
                    yp = empty_land(y(index), x(index), 4);
                    xp = empty_land(y(index), x(index), 3);
                    grid(yp, xp, 2) = 1;
                    grid(c.Y, c.X, 2) = 0;
                    c.X = xp;
                    c.Y = yp;
                    landed = true;
                else
                    % reverse course
                    c.Course = (5/6) * pi + c.Course + 1/3 * pi * rand;
                    c.DepartX = c.X;
                    c.DepartY = c.Y;
                end
            end
            
            if ~landed
                % chance to change course mid journey
                if rand < change_course_chance
                    c.Course = c.Course - pi / 2 + pi * rand;
                end
                
                num_empty_sea = sum(sum(empty_sea(:, :, 1)));
                
                if num_empty_sea > 0
                    pos = [];
                    diff = inf;
                    dist = 0;
                    backup_pos = [];
                    backup_pos2 = []; % this is a little ridiculous, I know
                    
                    for y = 1:3
                        for x = 1:3
                            if empty_sea(y, x, 1)
                                dx = empty_sea(y, x, 2) - c.DepartX;
                                dy = empty_sea(y, x, 3) - c.DepartY;
                                angle = atan2(dy, dx);
                                new_diff = abs(c.Course - angle);
                                new_dist = sqrt(dx ^ 2 + dy ^ 2);
                                backup_pos2 = empty_sea(y, x, 2:3);
                                
                                if new_dist > dist
                                    dist = new_dist;
                                    backup_pos = backup_pos2;
                                    
                                    if new_diff < diff && (dx ~= 0 || dy ~= 0)
                                        diff = new_diff;
                                        pos = backup_pos;
                                    end
                                end
                            end
                        end
                    end
                    
                    if isempty(pos)
                        pos = backup_pos;
                        
                        if isempty(pos)
                            pos = backup_pos2;
                        end
                    end
                    
                    c.PX = c.X;
                    c.PY = c.Y;
                    c.Times = 0;
                    grid(c.Y, c.X, 2) = 0;
                    grid(pos(2), pos(1), 2) = 1;
                    c.X = pos(1);
                    c.Y = pos(2);
                end
            end
            
            c.CanLand = max(c.CanLand - 1, 0);
        end

    	c = c.Next;
    end
    
    for j = 1:size(add_cells, 1)
        cells.add(add_cells(j, 1), add_cells(j, 2));
    end

    if mod(i - 1, render_interval) == 0
        if heat_output || bw_output
            orig_image(:, :, 1) = orig_image(:, :, 1) + uint8(grid(:, :, 2) * heat_value);
            orig_image(:, :, 2) = orig_image(:, :, 2) - uint8(grid(:, :, 2) * heat_value);
            orig_image(:, :, 3) = orig_image(:, :, 3) - uint8(grid(:, :, 2) * heat_value);
        end
        
        cell_mask = uint8(grid(:, :, 2));
        
        for j = 1:3
            if j == 1
                if ~bw_output; continue; end
                frame = orig_image;
            elseif j == 2
                if ~colour_output; continue; end
                frame = base_image;
            elseif j == 3
                frame = orig_image;
                if ~heat_output; continue; end
            end

            if j == 1
                frame(:, :, 1) = frame(:, :, 1) + cell_mask * 255;
                frame(:, :, 2) = frame(:, :, 2) + cell_mask * 255;
                frame(:, :, 3) = frame(:, :, 3) + cell_mask * 255;
            elseif j == 2
                uncell_mask = uint8(~cell_mask);
                frame(:, :, 1) = cell_mask * 255 + frame(:, :, 1) .* uncell_mask;
                frame(:, :, 2) = frame(:, :, 2) .* uncell_mask;
                frame(:, :, 3) = frame(:, :, 3) .* uncell_mask;
            end
            
            frame = insertText(frame, [5 5; 5 30], [i cells.Length]);
            
            if j == 1
                vid = video1;
            elseif j == 2
                vid = video2;
            elseif j == 3
                vid = video3;
            end
            
            writeVideo(vid, frame);
        end
        
        plot(population);
        drawnow;
        
        if pop_output
            writeVideo(video4, getframe);
        end
    end
    
%     if i > 1 && sound_output
%         pop_diff = population(i) - population(i - 1);
%         sound_indices = ((i - 1) * bits_per_step + 1):(i * bits_per_step);
%         sound_t = sound_indices / bits_per_step / steps_per_sec;
%         waveform = sin(2 * pi * sound_t * pop_diff * 150);
%         
%         if pop_diff < 0
%             waveform = sign(waveform); % square wave
%         end
%         
%         sound(sound_indices, 1) = waveform;
%     end
end

if bw_output
    close(video1);
end

if colour_output
    close(video2);
end

if heat_output
    close(video3);
end

if pop_output
    close(video4);
end

if sound_output
    for i = 1:steps
        if i > 1 && sound_output
            pop_diff = 0;
            avg_range = max(i - 20, 2):min(i + 20, steps);
            
            for j = avg_range
                pop_diff = pop_diff + abs(population(j) - population(j - 1));
            end
            
            pop_diff = pop_diff / length(avg_range);
            sound_indices = ((i - 1) * bits_per_step + 1):(i * bits_per_step);
            sound_t = sound_indices / bits_per_step / steps_per_sec;
            ang_freq = 2 * pi * sound_t * pop_diff * 200;
            waveform = sin(ang_freq) - floor(sin(ang_freq));
            waveform = waveform + (sin(ang_freq + 300) - floor(sin(ang_freq + 300))) * 0.3;
            if pop_diff < 0
                waveform = (waveform + sign(waveform)) / 2; % square wave
            end

            sound(sound_indices, 1) = waveform;
        end
    end

    soundsc(sound, sound_fs);
end

imshow(frame);
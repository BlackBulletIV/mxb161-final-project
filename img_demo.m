max_fertile_dist = 22; % max distance from sea before 0 fertility (change this for different images)
load_fertility = true; % whether to load fertility from file
map_noise = 10; % random noise in base map image

% display/output
green = [107 238 67]; % colour of fertile land
brown = [131 80 15]; % colour of barren land
sea = [25 95 220]; % colour of sea

%% Generate Grid

% load image
orig_image = imread('earth_img_demo.png');
figure
imshow(orig_image)
grid = double(orig_image);
grid_size = size(grid);

% this just uses red channel for now; need grayscale with typecasting out of
% uint8 to work for other images
grid(:, :, 1) = grid(:, :, 1) < 120; % land/sea
grid(:, :, 2) = 0; % cell present - 0 for now

half_height = grid_size(1) / 2;
sea_borders = [];

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

imshow(grid(:, :, 1) .* grid(:, :, 3))

y_mat = repmat((1:grid_size(1))', 1, grid_size(2));
x_mat = repmat(1:grid_size(2), grid_size(1), 1);

% temperature relative to equator (via y axis)
grid(:, :, 4) = max(min(1 - abs(half_height - y_mat)*1.15 / half_height, 1), 0);

imshow(grid(:, :, 4))

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

imshow(base_image)


% Create temperature lookup table (average temps in Celsius)
function temp = get_seasonal_temperature(location, season)
    % Temperature data structure
    temp_data = {
        % Location      Spring  Summer  Autumn  Winter
        'California',    17,     24,     17,     9;
        'Texas',         21,     29,     21,     12;
        'Minnesota',     8,     22,     9,     -7;
        'NewYork',       12,     24,     11,      2;
    };
    
    % Season mapping
    season_cols = {'spring', 'summer', 'autumn', 'winter'};
    
    % Handle autumn1 and autumn2 (both use autumn temps)
    if contains(season, 'autumn')
        season = 'autumn';
    end
    
    % Find indices
    loc_idx = find(strcmp(temp_data(:,1), location));
    season_idx = find(strcmp(season_cols, season)) + 1; % +1 for column offset
    
    temp = temp_data{loc_idx, season_idx};
end
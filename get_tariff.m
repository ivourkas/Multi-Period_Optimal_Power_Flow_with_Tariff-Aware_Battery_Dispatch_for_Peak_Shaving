function [C_sb, C_if] = get_tariff(tariffs, location, tariff_season, day_type, sector)
    try
        C_sb = tariffs.(location).(tariff_season).(day_type).(sector);
    catch
        error('Tariff not found for %s %s %s %s', location, tariff_season, day_type, sector);
    end
    C_if = C_sb;
end
function vrw = convertVRW(input,vrw)

switch input
    case 'ug/sqrt(Hz)'
        % Input is micro G per square root Hz
        % Output is m/s/sqrt(hr)

        vrw = vrw*1e-6*60*9.8;

    case 'm/s/sqrt(hr)'
        % Input is m/s/sqrt(hr)
        % Output is micro G per square root Hz
        % TIP: divide again by 60 to get m/s/sqrt(s)

        vrw = vrw/1e-6/60/9.8;

end

end
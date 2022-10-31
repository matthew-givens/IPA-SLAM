function arw = convertARW(input,arw)

switch input
    case 'rad/sqrt(s)'
        % Input is rad/sqrt(s)
        % Output is deg/sqrt(hr)

        arw = arw/pi*180*60;

    case 'deg/sqrt(hr)'
        % Input is deg/sqrt(hr)
        % Output is rad/sqrt(s)

        arw = arw*pi/180/60;

end

end
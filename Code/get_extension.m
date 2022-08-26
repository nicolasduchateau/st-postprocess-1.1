function extension = get_extension(number)

if ( (number >= 0) && (number < 10) )
    extension = ['_000',num2str(number)];
else
    if ((number >= 10) && (number < 100))
        extension = ['_00',num2str(number)];
    else
        if ((number >= 100) && (number < 1000))
            extension = ['_0',num2str(number)];
        else
            error('Too many frames inside Dicom series (> 1000 frames)');
        end
    end
end

end         % function get_extension

function dataPath = lojGetDataPath

    if ismac && strcmp(getenv('LOGNAME'), 'luis')
        dataPath = '/Users/luis/Box/boxPHD/Toolbox/2pdata/';
    elseif ispc && strcmp(getenv('username'), 'Luis')
        dataPath = 'C:\Users\Luis\Box Sync\boxEXPER\Toolbox\2pdata\';
    elseif ispc && strcmp(getenv('username'), 'dario')
        dataPath = 'C:\2pdata\';
    else
        error('unknown computer')
    end

end
%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% Type to iType
% 0:V,1:I,2:R,3:L,4:C
%--------------------------------------------------------------------------
function [Type] = funSimiType2Type(iType)
    switch iType
        case 0
            Type = 'V';
        case 1
            Type = 'I';
        case 2
            Type = 'R';
        case 3
            Type = 'L';
        case 4
            Type = 'C';
    end
end
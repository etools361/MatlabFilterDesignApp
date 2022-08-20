%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% Type to iType
% 0:V,1:I,2:R,3:L,4:C
%--------------------------------------------------------------------------
function [iType] = funSimType2iType(Type)
    switch Type
        case 'V'
            iType = 0;
        case 'I'
            iType = 1;
        case 'R'
            iType = 2;
        case 'L'
            iType = 3;
        case 'C'
            iType = 4;
    end
end
%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-19(yyyy-mm-dd)
% Suffix 2 double
%--------------------------------------------------------------------------
function DataOut=Suffix2Data(strDataIn)
    strUnit = {'y', 'e-24'; 'z', 'e-21';  'a', 'e-18'; 'F', 'e-15'; 'p',  'e-12'; 'N', 'e-9'; 'u',   'e-6'; 'k', 'e+3'; 'K', 'e+3';...
        'M', 'e+6'; 'Meg',   'e+6'; 'G',  'e+9'; 'T', 'e+12'; 'm', 'e-3';   'U',  'e-6'; 'n',  'e-9'; 'P', 'e+15'; 'f', 'e-15'; 'A', 'e-10'; ...
           'E', 'e+18'; 'Z', 'e+21'; 'Y', 'e+24'};
       
    N  = size(strUnit, 1);
     strDataIn = num2str(strDataIn);
     if strcmpi(strDataIn, 'inf')
         DataOut = inf;
         return;
     end
%      DataOut = str2double(strDataIn);
    flag = 0;
%     if isnan(DataOut)
        %% ??????????????
        for i=1:N
            strDataIn = regexprep(strDataIn, strUnit{i, 1}, strUnit{i, 2});
        end
        DataOut = str2num(strDataIn);

end


%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% [cellDataOut, Unit, strUnit0]=Data2Suffix(strDataIn);
%--------------------------------------------------------------------------
function [cellDataOut, Unit, strUnit0]=Data2Suffix(strDataIn, varargin)
    if nargin < 2
        pformat = '3.6';
    else
        pformat = varargin{1};
        if ischar(pformat)
        else
            if pformat == floor(pformat)
                pformat = sprintf('3.%d', pformat);
            else
                pformat = sprintf('%d.%d', floor(pformat), floor((pformat-floor(pformat))*10));
            end
        end
    end
    Unit = 1;
    strUnit0 = '';
    strUnit = {'y', 1e-24; 'z', 1e-21;  'a', 1e-18; 'F', 1e-15; 'p',  1e-12; 'N', 1e-9; 'u',   1e-6; 'k', 1e+3; 'K', 1e+3;...
        'M', 1e+6; 'Meg',   1e+6; 'G',  1e+9; 'T', 1e+12; 'm', 1e-3;   'U',  1e-6; 'n',  1e-9; 'P', 1e+15; 'f', 1e-15; 'A', 1e-10; ...
           'E', 1e+18; 'Z', 1e+21; 'Y', 1e+24};
    strUnit2 = {'y', 1e-24; 'z', 1e-21; 'a', 1e-18; 'f', 1e-15; 'p', 1e-12; 'n',  1e-9; 'u',  1e-6; 'm',  1e-3; '', 1;...
                'K',  1e+3; 'M', 1e+6;  'G',  1e+9; 'T', 1e+12; 'P', 1e+15; 'E', 1e+18; 'Z', 1e+21; 'Y', 1e+24;};
       
    N  = size(strUnit, 1);
    N2 = size(strUnit2, 1);
    mData = length(strDataIn);
    for ii=1:mData
        strDataIn0 = strDataIn(ii);
        try
            strDataIn0 = num2str(strDataIn0);
            DataOut = str2double(strDataIn0);
            flag = 0;
            if isnan(DataOut)
                %% ?????????????????????
                for i=1:N
                    if strfind(strDataIn0, strUnit{i, 1})
                        strTemp = strDataIn0;
                        strTemp   = regexp(strTemp, strUnit{i, 1}, 'split');
                        DataOut = str2double(strTemp{1})*strUnit{i, 2};
                        if isnan(DataOut)
                            warning('0x:%s format error!\n', strDataIn0);
                            strDataOut = strDataIn0;
                        else
                            flag = 1;
                        end
                        break;
                    end
                end
                if flag == 0
                    warning('1x:%s format error!\n', strDataIn0);
                    strDataOut = strDataIn0;
                else
                end
            end

            %% ???????????????
            PN = 1;
            if DataOut == 0
                PN = 0;
                strDataOut = '0';
            elseif DataOut <0
                PN = -1;
                DataOut = -DataOut;
            elseif DataOut == Inf
                PN = 0;
                strDataOut = 'Inf';
            end
            if PN ~= 0
                for i=1:N2
                    if DataOut < strUnit2{i, 2}
                        if i<2
                            i = 2;
                        end
                        Unit    = strUnit2{i-1, 2};
                        strUnit0 = strUnit2{i-1, 1};
                        Datax = DataOut/strUnit2{i-1, 2};
                        if nargin < 2
                            strNum = char(vpa(Datax,100));
                            n = length(strNum);
                            m = n-strfind(strNum, '.');
                            strp = sprintf('%%%d.%df', n-m-1, m);
                        else
                            strp = sprintf('%%%sf', pformat);
                        end
                        strSt = sprintf(strp, Datax);
                        NstrSt = length(strSt);
                        if ~isempty(strfind(strSt, '.'))
                            for j = 1:NstrSt
                                if strcmp(strSt(end-j+1), '0')
                                    strSt2 = '0';
                                else
                                    if strcmp(strSt(end-j+1), '.')
                                        strSt2 = strSt(1:end-j+1-1);
                                    else
                                        strSt2 = strSt(1:end-j+1);
                                    end
                                    break;
                                end
                            end
                        else
                            warning('2x:%s format error!\n', strSt);
                            strSt2 = strSt;
                        end
                        strDataOut = sprintf('%s%s', strSt2, strUnit2{i-1, 1});
                        break;
                    else
                    end
                end
            end
            if PN < 0
                strDataOut = ['-' strDataOut];
            end
        catch
            if iscell(strDataIn0)
                strDataOut = char(strDataIn0);
            else
                strDataOut = strDataIn0;
            end
        end
        cellDataOut0{ii} = strDataOut;
    end
    if ii==1
        cellDataOut = cellDataOut0{1};
    else
        cellDataOut = cellDataOut0;
    end
end


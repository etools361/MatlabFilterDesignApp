%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 找出带节点的桥
%--------------------------------------------------------------------------
function [iBn, cellLoop, cellDevice, iLoop, Looping, e] = funFindBridgeWiNode(h, e, Bn, iBn)
% 排除GND
Bn(Bn==1)  = [];
cellLoop   = [];
cellDevice = [];
iLoop      = 0;
Looping    = 0;
% 找出剩余节点，这些节点都是桥，首先找出节点桥
for ii=Bn
    [iBn, cellLoop, cellDevice, iLoop, Looping, e] = funDFSFindBridge(h, e, ii, iBn, cellLoop, cellDevice, iLoop, Looping);
end
end
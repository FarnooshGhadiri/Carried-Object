function [winMtx_mask, winThreshold,bias] = winnerPopOut(VR,nonZeros,method)
%[winMtx_mask, winThreshold,bias] = winnerPopOut(VR,nonZeros,method,bias)
% choose the most potential winner
%
%
% Liming Wang, Jan 2008
%


% mn is the mean value of 'VR'
% 1) nonZeros==0
%       only computer 'mn' on the region of 'VR' where VR>0
% 2) nonZeros!=0
%       computer 'mn' all over 'VR'

if(~exist('nonZeros','var'))
    nonZeros    = 1;
end
if(~exist('method','var'))
    method=2;
end

if(nonZeros==0)
    mn = sum(sum(VR))/prod(size(VR));
else
    mn = sum(sum(VR))/(sum(sum(VR>0))+eps);
end

winMtx_mask = [];
winThreshold = 10;
if (method ==1)
% use first moment
% those point which is greater than (bias+1)*mn is selected
% here we computer bias according the max(VR)
    m_bias = max(max(VR))/mn;
    % here, 0.5 is just a imperical value,I have not prove its correctness
    bias = 0.5*m_bias;
    winThreshold = (bias+1)*mn;
    winMtx_mask = VR>winThreshold;    
else
    if (method ==2)
        % use second moment
        stdd=stdDev(VR,mn,nonZeros);
        m_bias = (max(max(VR)) - mn)/(stdd+eps);
        bias = 0.4*m_bias; 
        winThreshold = mn + bias*stdd;
        winMtx_mask = VR>winThreshold;        
    elseif(method == 3)
        % use zero moment
        m_bias = max(max(VR));
        bias = 0.5*m_bias;
        winMtx_mask = VR>bias;
        winThreshold = bias;
    end
end

function res = stdDev(mtx,mn,nonZeros)
if(nonZeros==0)
    res = sum(sum(abs(mtx-mn).^2)) / max((prod(size(mtx)) - 1),1);
else
    res = sum(sum((abs(mtx-mn).^2).*(mtx>0))) / max((sum(sum(mtx>0)) - 1),1);
end
res = sqrt(res);

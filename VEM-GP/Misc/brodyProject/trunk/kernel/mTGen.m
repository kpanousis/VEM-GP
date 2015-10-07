function [m_t,varargout] = mTGen(times,clickTimes)
%MTGEN counts the number of clicks which have occurred at any time.
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%28 February, 2014
%mTGen takes two sorted lists, times and clickTimes, and returns m_t, the
%number of entries in clickTimes which are smaller than each entry in
%times.  May also return tfh, the time each click was first heard, i.e.,
%entry i in tfh is  the first index in times which corresponds to an entry
%which is strictly greater than clickTimes(i)
nExtraOutputs = max(0,nargout - 1);
if numel(clickTimes) == 0;
    m_t = zeros(size(times));
    tfh = [];  %Is this the right choice?
elseif times(end) < clickTimes(1)
    m_t = zeros(size(times));
    if nExtraOutputs > 0
        tfh = nan(size(clickTimes));
    end
elseif times(1) > clickTimes(end)
    m_t = numel(clickTimes) * ones(size(times));
    if nExtraOutputs > 0
        tfh = ones(size(clickTimes));
    end
else %the lists overlap
    %Preallocate m_t
    m_t = zeros(size(times));
    if nExtraOutputs > 0
        tfh = nan(size(clickTimes));
    end
    clicksSoFar = 0;
    for timeIdx = 1:numel(times)
        if clicksSoFar == numel(clickTimes)
            m_t(timeIdx:end) = clicksSoFar;
            break
        end
        while numel(clickTimes) > clicksSoFar && times(timeIdx) > clickTimes(clicksSoFar + 1)
            clicksSoFar = clicksSoFar + 1;
            if nExtraOutputs > 0
                tfh(clicksSoFar) = timeIdx;
            end
        end
        m_t(timeIdx) = clicksSoFar;
    end
end


if nExtraOutputs > 0
    varargout{1} = tfh;
end
end
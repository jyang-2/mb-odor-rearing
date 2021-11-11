function [riseTime, fallTime, midTime] = get_frameout_times(Frame_Out, time)
    frameOutLogical = logical(Frame_Out);
    frameOutDiff = diff(frameOutLogical);
    risingEdge = frameOutDiff>0;
    fallingEdge = frameOutDiff<0;

    riseTime = time(risingEdge);
    fallTime = time(fallingEdge);
    midTime = mean([riseTime fallTime],2)'; 
end
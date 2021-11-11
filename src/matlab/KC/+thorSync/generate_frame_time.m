function [frameTime, indexes] = generate_frame_time(Frame_Out, time)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%% below is the new vision
frameOutLogical = logical(Frame_Out);
frameOutDiff = diff(frameOutLogical);
risingEdge = find(frameOutDiff>0);
fallingEdge = find(frameOutDiff<0);
len =fallingEdge - risingEdge;
maxLen = max(len);
minLen = min(len);
frameOutDiff = diff(Frame_Out);
% if maxLen > 1.5*minLen
%     threshold = minLen + (maxLen - minLen)/2;
%     frameOutDiff(risingEdge(len>threshold))=0;
% end
frameOutDiff = vertcat(0,frameOutDiff);

z1 = frameOutDiff;
indexes = find(z1>0);
frameTime = time(indexes);



end


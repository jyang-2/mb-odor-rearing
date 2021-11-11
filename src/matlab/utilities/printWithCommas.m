function str = printWithCommas(a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
str = [sprintf('%s,', a(1:end-1)), sprintf('%s', a(end))];
end


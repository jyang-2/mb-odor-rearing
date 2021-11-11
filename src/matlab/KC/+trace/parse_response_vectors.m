function [response_ind, response_peak] = parse_response_vectors(bin_response,peak_ind, peak_response)
%PARSE_RESPONSE_VECTORS Summary of this function goes here
%   Detailed explanation goes here

response_ind = nan(size(bin_response));
response_peak = nan(size(bin_response));

response_ind(bin_response) = peak_ind(bin_response);
response_peak(bin_response) = peak_response(bin_response);


end


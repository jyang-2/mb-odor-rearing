function myuuid = generate_uuid()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

temp =  java.util.UUID.randomUUID;
myuuid = temp.toString;


end


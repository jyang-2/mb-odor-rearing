function ver = get_xml_software_version(exml)
%GET_XML_SCOPETYPE Summary of this function goes here
%   Detailed explanation goes here
ver = str2double(strsplit(exml.Software.versionAttribute, '.'));
end


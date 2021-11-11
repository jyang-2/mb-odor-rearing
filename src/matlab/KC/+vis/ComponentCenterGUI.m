classdef ComponentCenterGUI < matlab.graphics.chartcontainer.ChartContainer
    %COMPONENTCENTERGUI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cm    (:,3) double = Nan
        siz   (1,3)
        XData (:,1) double = NaN
        YData (:,1) double = NaN
        ZData (:,1) double = NaN
    end
    
    methods
        function obj = ComponentCenterGUI(inputArg1,inputArg2)
            %COMPONENTCENTERGUI Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


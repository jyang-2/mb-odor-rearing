for z = 1:d3
    tb = axtoolbar(harray(z), {'export', 'pan', 'zoomin', 'zoomout', 'restoreview'});
    btn = axtoolbarbtn(tb, 'push');
    btn.Icon = 'roi_icon.png';
    btn.Tooltip = 'draw ROI';
    btn.ButtonPushedFcn = @customcallback; 
end

function customcallback(src,event)
        drawfreehand(event.Axes);
    end
function fig = plot_thorsync_lines(ai)
    if isfield(ai, 'FrameOut')
        ai.Frame_Out = ai.FrameOut;
        ai = rmfield(ai, 'FrameOut');
    end
    fig = figure;
    plot(ai.time, ai.scopePin); hold on; 
    plot(ai.time, ai.olfDispPin);
    plot(ai.time, ai.Frame_Out);
    xlabel('time (sec)');
    ylabel('V');
    ylim([-1 6])
    title('thorsync lines', 'Interpreter', 'none');
    legend({'scopePin', 'olfDispPin', 'Frame_Out'}, 'Location', 'best', 'Interpreter', 'none');
    legend('boxoff');
end
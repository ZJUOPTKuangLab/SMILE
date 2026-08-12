function slider_and_button_demo1(tmp_data, angle_list)
    % create figure and subplot
    fig = figure;
    ax = subplot(1,1,1);
    num = numel(angle_list);
    delta = angle_list(2) - angle_list(1);
    stepsize = 1/(num-1);
    

    temp = tmp_data{1, (num+1)/2};
    tmp_x = temp(:,1);
    N11 = temp(:,2);
    % create initial plot
    scatter(tmp_x, N11,1.5,'filled');
    xlim([min(tmp_x)-10, max(tmp_x)+10])
    title(angle_list((num+1)/2))
    
    % create slider
    slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', num, 'Value', (num+1)/2, ...
                       'Position', [100 20 120 20],'SliderStep', [stepsize, stepsize*5], 'Callback', @slider_callback);
    
    
    % slider callback function
    function slider_callback(source, event)
        % get slider value
        slider_value = get(source, 'Value');
        
        % update plot data
        temp1 = tmp_data{1, slider_value};
        tmp_x1 = temp1(:,1);
        N111 = temp1(:,2);
        scatter(tmp_x1, N111,1.5,'filled');
        xlim([min(tmp_x)-10, max(tmp_x)+10])
        title(angle_list(slider_value))
        


    end

end
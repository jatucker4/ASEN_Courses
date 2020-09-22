function [] = Plot_Settings(psi_plot,phi_plot,press_plot,c,alpha_vec,V_inf_vec,flag)
%Plot_Settings Performs handling of the figure objects and creates the plot
%settings necessary for the data to be displayed in a concise manner.
%
% Author: Johnathan Tucker
% Collaborators: N/A
% This function takes in the figure objects as well, the altered variables
% detailed in the third bullet point, and a switch flag
%
% Last Revised: 2/27/2020
%% Begin Plot Setting
% Switch through each case. This first block is for the first bullet point
if flag == 1
    % Create figure
    figure(14);
    % Create first subplot with plot settings
    s1 = subplot(1,2,1);
    title('$Stream\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 4 -2 2])
    colorbar
    % Create second subplot with plot settings
    s2 = subplot(1,2,2);
    title('$Potential\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 4 -2 2])
    colorbar
    sgtitle('$Baseline\:Potential\:and\:Stream\:Functions$','Interpreter','latex')
    % Access the data children of the figure objects
    psi_child = get(psi_plot.Children,'children');
    phi_child = get(phi_plot.Children,'children');
    % Determine the data type of the child and handle it accordingly
    % copyobj takes the data from the child and puts it into the subplot
    if iscell(phi_child)
        copyobj(phi_child{end},s2);
    else
        copyobj(phi_child,s2);
    end
    if iscell(psi_child)
        copyobj(psi_child{end},s1);
    else
        copyobj(psi_child,s1);
    end
    % Maximize
    set(gcf, 'Position', get(0, 'Screensize'));
    % Create pressure plot for the first bullet points
    figure(15)
    ax1 = gca;
    title('$Baseline\:Pressure\:Plot$','Interpreter','latex','FontSize',23)
    xlabel('$x-distance\:[m]$','Interpreter','latex','FontSize',13)
    ylabel('$y-distance\:[m]$','Interpreter','latex','FontSize',13)
    axis([-2 4 -2 2])
    colorbar
    press_child = get(press_plot.Children,'children');
    if iscell(press_child)
        copyobj(press_child{end},ax1);
    else
        copyobj(press_child,ax1);
    end
    % Maximize
    set(gcf, 'Position', get(0, 'Screensize'));
    
% Begin block for the chord sensitivity test
elseif flag == 2 || flag == 3 || flag == 4 || flag == 5 || flag == 6 || flag == 7
    figure(15+flag-1)
    % Create first subplot with plot settings
    s1 = subplot(1,2,1);
    title('$Stream\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 flag+2 -2 2])
    colorbar
    % Create second subplot with plot settings
    s2 = subplot(1,2,2);
    title('$Potential\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 flag+2 -2 2])
    colorbar
    sgtitle(sprintf('$Stream\\:and\\:Potential\\:Contours\\:With\\:c\\:=\\:%d\\:m$',c),'Interpreter','latex')
    % Access the data children of the figure objects
    psi_child = get(psi_plot.Children,'children');
    phi_child = get(phi_plot.Children,'children');
    % Determine the data type of the child and handle it accordingly
    % copyobj takes the data from the child and puts it into the subplot
    if iscell(psi_child)
        copyobj(psi_child{end},s1);
    else
        copyobj(psi_child,s1);
    end
    if iscell(phi_child)
        copyobj(phi_child{end},s2);
    else
        copyobj(phi_child,s2);
    end
    % Maximize
    set(gcf, 'Position', get(0, 'Screensize'));
    
% Begin block for the angle of attack sensitivity test
elseif flag == 8 ||flag == 9||flag == 10||flag == 11||flag == 12||flag == 13
    figure(21+flag-7)
    % Create first subplot with plot settings
    s1 = subplot(1,2,1);
    title('$Stream\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 4 -2 2])
    colorbar
    % Create second subplot with plot settings
    s2 = subplot(1,2,2);
    title('$Potential\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 4 -2 2])
    colorbar
    sgtitle(sprintf('$Stream\\:and\\:Potential\\:Contours\\:With\\:\\alpha\\:=\\:%d\\:Degrees$',alpha_vec),'Interpreter','latex')
    % Access the data children of the figure objects
    psi_child = get(psi_plot.Children,'children');
    phi_child = get(phi_plot.Children,'children');
    % Determine the data type of the child and handle it accordingly
    % copyobj takes the data from the child and puts it into the subplot
    if iscell(phi_child)
        copyobj(phi_child{end},s2);
    else
        copyobj(phi_child,s2);
    end
    if iscell(psi_child)
        copyobj(psi_child{end},s1);
    else
        copyobj(psi_child,s1);
    end
    % Maximize
    set(gcf, 'Position', get(0, 'Screensize'));
    
% Begin block for the free-stream velocity sensitivity test
elseif flag == 14 ||flag == 15||flag == 16||flag == 17||flag == 18||flag == 19
    figure(22+flag)
    % Create first subplot with plot settings
    s1 = subplot(1,2,1);
    title('$Stream\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 4 -2 2])
    colorbar
    % Create second subplot with plot settings
    s2 = subplot(1,2,2);
    title('$Potential\:Function$','Interpreter','latex')
    xlabel('$x-distance\:[m]$','Interpreter','latex')
    ylabel('$y-distance\:[m]$','Interpreter','latex')
    axis([-2 4 -2 2])
    colorbar
    sgtitle(sprintf('$Stream\\:and\\:Potential\\:Contours\\:With\\:V_{inf}\\:=\\:%d\\:\\frac{m}{s}$',V_inf_vec),'Interpreter','latex')
    % Access the data children of the figure objects
    psi_child = get(psi_plot.Children,'children');
    phi_child = get(phi_plot.Children,'children');
    % Determine the data type of the child and handle it accordingly
    % copyobj takes the data from the child and puts it into the subplot
    if iscell(phi_child)
        copyobj(phi_child{end},s2);
    else
        copyobj(phi_child,s2);
    end
    if iscell(psi_child)
        copyobj(psi_child{end},s1);
    else
        copyobj(psi_child,s1);
    end
    % Maximize
    set(gcf, 'Position', get(0, 'Screensize')); 
end

end

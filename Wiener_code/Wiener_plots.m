%% plot for the Wiener Filter code 
% plots for the Wiener part of the papaer, so all the predicitons, and then
% only flexor and extensor
% saving it somewhere 

% complete plot 
figure
time = linspace(1,(size(Y_fin,1))*1/new_samp,size(Y_fin,1));
for k = 1:size(Y_test,2)   
    sub(k) = subplot(size(Y_test,2),1,k);
    plot(time,Y_test(:,k),'b','LineWidth',1.1)
    hold on
    plot(time,Y_fin(:,k),'r','LineWidth',1.1)
    xlim([20 40])
    set(gca,'Fontsize',18);
end


% only the two muscles (flexor and extensor)
h = figure;
to_plot = [5];
to_plot_name = {'example 1 '};
for i = 1:length(to_plot)
    sub(i) = subplot(length(to_plot),1,i);
    plot(time,Y_test(:,to_plot(i)),'b','LineWidth',1.1)
    ylabel(to_plot_name{i})
    hold on
    plot(time,Y_fin(:,to_plot(i)),'r','LineWidth',1.1)
    box('off')
    set(gca,'FontSize',18)
end
linkaxes(sub,'x')


%

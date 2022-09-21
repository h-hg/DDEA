function plot_position
t1=plot(archive_save_pos(8,:,1),archive_save_fit(8,:,1),'*');
set(t1,'LineWidth',2);
hold on
t2=plot(current_save_pos(8,:,1),current_save_fit(8,:,1),'rp');
set(t2,'LineWidth',2);
hold on
hx1=xlabel('Position on the 1-th dimension');
set(hx1,'FontSize',16);
hx2=ylabel('Fitness value');
set(hx2,'FontSize',16);
set(gca,'FontSize',10);
end
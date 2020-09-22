function plot_truss(nodetable, connecttable)

figure
hold on

for m=1:size(connecttable,1)
    plot3([nodetable(connecttable(m,2),2);nodetable(connecttable(m,3),2)],...
          [nodetable(connecttable(m,2),3);nodetable(connecttable(m,3),3)],...
          [nodetable(connecttable(m,2),4);nodetable(connecttable(m,3),4)],'k');
end

hold off
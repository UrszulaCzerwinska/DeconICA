function plotAverageStability(folder)

ls = dir(folder);

count = 1;

for i=3:size(ls,1)
    name = ls(i).name;
    k = size(strfind(name,'stability.txt'),1);
    if(k>0)
        xx = load(strcat(folder,name));
        array(count,1) = size(xx,1);
        array(count,2) = mean(xx);
        count = count+1;
    end
end

fullFileName2 = fullfile(folder, '/avg_stability.plot.txt');


T=[array(:,1),array(:,2)];

dlmwrite(fullFileName2 ,T, '\t');

fig_stab = bar(array(:,1),array(:,2)); ylim([0 1]);


fullFileName = fullfile(folder, '/stability.png');

saveas(fig_stab,fullFileName);


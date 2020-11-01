function [outputArg1,outputArg2] = phase_space(t,dt,r,v)
 if mod(ceil(t/dt),ceil(5/dt)) == 0
        figure
        plot(r,v,'o','MarkerSize', 2)
        xlabel('X')
        ylabel('Y')
        tit = 'Electron Phase-space distribution for t = ';
        tit = strcat(tit,num2str(t));
        title(tit);
        fname = strcat(strcat('phase_space_',num2str(t)),'.png')
        saveas(gcf,fname)
 end
end


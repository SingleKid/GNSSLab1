
step =4;
data = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\out1.txt');
refblh = [51.258643*pi/180     -114.100492*pi/180  	1127.345];
blh = xyz2blh(data);
blh(:,1) = blh(:,1) * 180 / pi;
blh(:,2) = blh(:,2) * 180 / pi;

if step == 1
    l1 = scatter(blh(:,2), blh(:,1), '.');
    hold on;
    l2 = scatter(refblh(:,2), refblh(:,1), 'o', 'r', 'linewidth',4);
    legend([l1, l2], {'Samples', 'Truth'});
    grid on;
    box on;
    xlabel('Longitude / °');
    ylabel('Latitude / °');
    set(gca,'FontSize',16);
end


% mean_blh = mean([refblh; blh(:,:)]);
enu = xyz2enu(refblh, data);

enu_error = enu;
std_ee = std(enu_error)';
rms_ee = rms(enu_error)';
mean_ee = mean(enu_error)';

if step == 2
    l1 = plot(blh(:,3));
    hold on;
    l2 = plot([-500, 4000], [refblh(:,3), refblh(:,3)]);
    legend([l1, l2], {'Samples', 'Truth'});
    grid on;
    box on;
    xlabel('Epoch time / seconds');
    ylabel('Altitude / meters');
    set(gca,'FontSize',16);
end

if step == 3
    env = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\env.txt');
    yl = {'East Error / meters'; 'North Error / meters'; 'Up Error / meters'};
%     range = [2.5,2,6];
    for i = 1 : 3
        subplot(3,1,i);
        hold on;
        scatter(1:3601,enu_error(:,i),'.');
        grid on;
        b = mean(enu_error(:,i));
        s = std(enu_error(:,i));
        axis([0, 3600, min(enu_error(:,i) - env(:,i))-0.2, max(enu_error(:,i) + env(:,i))+0.2]);
        
        plot(1:3601,enu_error(:,i) + env(:,i),'r');
        plot(1:3601,enu_error(:,i) - env(:,i),'r');

%         scatter(1:3601,enu_error(:,i) + env(:,i),'r','.');
%         scatter(1:3601,enu_error(:,i) - env(:,i),'r','.');
%         plot([0, 3600], [b + s, b + s], 'r');
%         plot([0, 3600], [b - s, b - s], 'r');
        ylabel(yl{i});
        box on;
    end
    xlabel('Epoch time / seconds');
end

if step == 4
    hold on;
    dops = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\dops.txt');
    a = [plot(dops(:,1)), plot(dops(:,2)), plot(dops(:,3)), plot(dops(:,4)), plot(dops(:,5))];
    legend(a, {'EDOP', 'NDOP', 'VDOP', 'HDOP', 'PDOP'});
    axis([-400, 4000, 0, max(dops(:,5))+1]);
    grid on;
    box on;
    xlabel('Epoch time / seconds');
    ylabel('DOP Value / meters');
    set(gca,'FontSize',16);
end

if step == 5
    sat_num = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\satn.txt');
    plot(sat_num);
    grid on;
    box on;
    xlabel('Epoch time / seconds');
    ylabel('Number of Staellites');
    mins = min(sat_num);
    maxs = max(sat_num);
    axis([-400, 4000, mins-1, maxs+1]);
end

if step == 6
%     Get the elevation-residual relationship
    el = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\res2.txt');
    el = el .* 180 / pi;
    re = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\log.txt');
    hold on;
    
    sat_num = 0;
    legends = {''};
    handles = [];
    for i = 1 : 32
        if(range(re(:,i)) ~= 0)
            s = scatter(el(:,i), re(:,i),300, '.');
            sat_num = sat_num + 1;
            handles(sat_num) = s;
            legends{sat_num} = sprintf('PRN %2d', i);
        end
    end
    
    
    for i = 1 : sat_num
        set(handles(i), 'MarkerEdgeColor', [rand, rand, rand]);
    end
    legend(handles, legends);
%     e = abs(data(:,1) * 180 / pi);
%     
%     r = data(:,2);
% 
%     scatter(e, r, '.');
    grid on;
    box on;
    xlabel('Elevation / degrees');
    ylabel('Residual / meters');
    axis([0,90, -5,3])
    set(gca,'FontSize',16);
end

if step == 7
    res = load('C:\Users\80575\OneDrive\文档\Visual Studio 2015\Projects\GNSSLab1\GNSSLab1\log.txt');
    nepoch = size(res,1);
    sat_num = 0;
    epoch_num = size(res,1);
    prn_list = zeros(32, 0);
    res_t = zeros(epoch_num, 32);
    for i = 1 : 32
        if range(res(2:end,i)) ~= 0
            sat_num = sat_num + 1;
            res_t(:,sat_num) = res(:,i);
            prn_list(sat_num) = i;
        end
    end
    for i = 1 : sat_num
        subplot(sat_num, 1, i);
        hold on;
    
        last = 1;
        for j = 2 : epoch_num
            if res_t(j, i) == 0
                if res_t(last ,i) ~= 0
                    plot(last:j-1,res_t(last:j-1, i), 'r');
                end
                last = j;
            elseif res_t(last ,i) == 0
            last = last + 1;
            end
        end
    %     if last == 1
        plot(last:epoch_num,res_t(last:epoch_num,i), 'r');
        ylabel(sprintf('PRN %d', prn_list(i)));
    %     end
%         axis([0,nepoch,-1.5,1.5]);
        box on;
        grid on;
%      plot(res_t(:,i));
    end
end
% scatter3(enu(:,1), enu(:,2), enu(:,3), '.');
% kmlwrite('V:\a.kml', blh(:,1), blh(:,2))
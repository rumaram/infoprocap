classdef Plotter
    methods(Static)

        function Cm=capMat(obj,C,filename)
            arguments
                obj 
                C 
                filename ="no_save"
            end
            K=obj.K;
            %=====Calculating=====================
            C_tot=sum(C,"all");
            Cm=zeros(obj.max_deg+1);
            if size(obj.u,2)~=2
                disp("Error: Feature size should be 2 to generate Capacity matrix");
                return
            end

            for i=1:obj.basis_size               
                Cm(obj.degrees(i,1)+1,obj.degrees(i,2)+1)=C(i); 
            end
            %=====Plotting===========================
            if filename=="no_plot"
                return;
            elseif filename=="no_save"
                figure('Visible','on');
            else
                figure('Visible','off');
            end

            imagesc(Cm);

            n_color=50;
            color1   = [1, 1, 1];
            color2 = [0.0 0.5 0.0];
            customMap = [linspace(color1(1), color2(1), n_color)', ...
                linspace(color1(2), color2(2), n_color)', ...
                linspace(color1(3), color2(3), n_color)'];
            colormap(customMap);
            colorbar;

            [rows, cols] = size(Cm);

            border_color=[0.5,0.5,0.5];
            % Draw vertical lines
            for c = 0.5:1:cols+0.5
                line([c c], [0.5 rows+0.5], 'Color', border_color, 'LineWidth', 0.5);
            end

            % Draw horizontal lines
            for r = 0.5:1:rows+0.5
                line([0.5 cols+0.5], [r r], 'Color', border_color, 'LineWidth', 0.5);
            end

            % Overlay gray cells for values above diagonal
            mask = triu(true(size(Cm)), 1);
            mask=flipud(mask);
            [row, col] = find(mask);
            for i = 1:length(row)
                % Draw a black rectangle on top of the cell
                rectangle('Position', [col(i)-0.5, row(i)-0.5, 1, 1], ...
                    'FaceColor', border_color, 'EdgeColor', border_color);
            end

            %====Circle fraction===============

            pos = [0.63 0.7 0.14 0.14];   % tweak to taste (top-left)
            % Draw circle
            annotation('ellipse', pos, ...
                'Units','normalized', ...
                'FaceColor','w', ...
                'Color','k', ...
                'LineWidth',1);   % circle outline

            % Put fraction text centered inside circle
            str = sprintf('$\\frac{%g}{%g}$', round(C_tot,2), obj.K+1);
            annotation('textbox', pos, ...
                'Units','normalized', ...
                'String', str, ...
                'Interpreter','latex', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'EdgeColor','none', ...
                'FontSize',16);
            %==============================

            xlabel("Degree of u_1",'Interpreter', 'tex');
            ylabel("Degree of u_2",'Interpreter', 'tex');

            set(gca, 'YDir', 'normal')
            ax = gca;
            ax.FontSize = 14;
            set(gca, 'XTick', 1:size(Cm,2), 'XTickLabel', 0:size(Cm,2)-1);
            set(gca, 'YTick', 1:size(Cm,1), 'YTickLabel', 0:size(Cm,1)-1);
            ax.XLabel.FontSize=16;
            ax.YLabel.FontSize=16;
            axis square;

            set(gcf, 'Units', 'pixels', 'Position', [100 100 600 600]);

            if filename~="no_save"
                exportgraphics(gcf, filename, 'Resolution', 300);
            end

        end

        function Cd=capBar(obj,C,filename,y_lim)
            arguments
                obj 
                C 
                filename="no_save"
                y_lim =obj.K % y axis limit
            end
            % ===============Calculating===================
            C_tot=sum(C,"all");
            Cd=zeros(obj.max_deg+1,4);    
            % 1st column-> 1 element terms,
            % 2nd column-> 2 element terms
            % 3rd column-> more than 2 element terms
            % 4th column-> max total capacity

            for i=1:obj.basis_size
                current_degrees=obj.degrees(i,:);
                tot_degree=sum(current_degrees,"all");
                num_nonzeros=nnz(current_degrees);

                if num_nonzeros==0
                    Cd(1,1)=C(1);
                elseif num_nonzeros==1
                    Cd(tot_degree+1,1)=Cd(tot_degree+1,1)+C(i);
                elseif num_nonzeros==2
                    Cd(tot_degree+1,2)=Cd(tot_degree+1,2)+C(i);
                else
                    Cd(tot_degree+1,3)=Cd(tot_degree+1,3)+C(i);
                end
                Cd(tot_degree+1,4)=Cd(tot_degree+1,4)+1;

            end

            %============Plotting=======================
            if filename=="no_plot"
                return;
            elseif filename=="no_save"
                figure('Visible','on');
            else
                figure('Visible','off');
            end

            h=bar(Cd(:,1:3),'stacked');
            xlim([0.5, obj.max_deg+1 + 0.5]);        
            xticks(1:obj.max_deg+1);               % Bar positions
            xticklabels(0:obj.max_deg);        % Labels you want        
            ylim([0,y_lim]);

            h(1).FaceColor=[248, 237, 140]/255;
            h(2).FaceColor=[211, 230, 113]/255;
            h(3).FaceColor=[137, 172, 70]/255;
            legend(h, {'Single term', '2 terms', 'More than 2 terms'},'Location', 'northeast',FontSize=14);

            %==numbers above bars=========
            bar_tops = sum(Cd(:,1:3), 2);
            for i = 1:obj.max_deg+1
                Csum=round(sum(Cd(i,1:3)),2);
                bar_x = i;               % x-position of the bar
                bar_y = bar_tops(i);     % y-position (top of the bar)
                text(bar_x, bar_y+1, num2str(Csum), 'HorizontalAlignment', 'center',FontSize=14);
            end
            %====Circle fraction===============
            pos = [0.15 0.75 0.14 0.14];   % tweak to taste (top-left)

            % Draw circle
            annotation('ellipse', pos, ...
                'Units','normalized', ...
                'FaceColor','w', ...
                'Color','k', ...
                'LineWidth',1);   % circle outline

            % Put fraction text centered inside circle
            str = sprintf('$\\frac{%g}{%g}$', round(C_tot,2), obj.K+1);
            annotation('textbox', pos, ...
                'Units','normalized', ...
                'String', str, ...
                'Interpreter','latex', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'EdgeColor','none', ...
                'FontSize',18);
            %==============================
            xlabel("Degree of basis");
            ylabel("Capacity");
            axis square;
            ax = gca;
            ax.FontSize = 14; 
            ax.XLabel.FontSize=18;
            ax.YLabel.FontSize=18;
            set(gcf,'Units','pixels','Position',[100 100 600 600]);  % [left bottom width height]
            if filename~="no_save"
                exportgraphics(gcf, filename, 'Resolution', 300);
            end
        end

    end
end
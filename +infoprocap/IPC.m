classdef IPC<handle
    properties
        progress_display=true;
        u;          %inputs.
        u_basis;    % cell array of first max_deg no. of individual legendre basis functions (not products) calculated at input u
        y;          %calculated product basis functions of inputs u.        
        max_deg;
        degrees;
        basis_size;
        sample_size;
        dimn;
        threshs;
        samps_arr;
    end

    methods
        function obj=IPC(u,max_deg)
            obj.max_deg=max_deg;
            obj.u=u;
            obj.sample_size=size(u,1);
            obj.dimn=size(u,2);

            obj.degrees=IPC.stars_and_bars(obj.max_deg,obj.dimn);  % dsitribution of basis is equivalent to a stars and bars problem
            
            obj.basis_size=size(obj.degrees,1);
            obj.y=ones(obj.sample_size,obj.basis_size);

            obj.orthoPoly3();
            for basis=1:obj.basis_size
                for q=1:obj.dimn
                    deg=obj.degrees(basis,q);
                    obj.y(:,basis)=obj.y(:,basis).*obj.u_basis(:,q,deg+1);
                end             
            end

            if obj.progress_display
                disp("Basis Initiated");
            end

        end

        function initThresholds(obj,X)
            X=cat(2,X,ones(size(X,1),1));
            N=size(X,1);

            [Q,~]=qr(X,0);
            Xt=sqrt(N)*Q;
            
            term1=mean((sum(Xt.^2,2)).^2,1);  
            for bsis=1:obj.basis_size
                
                term2=1;
                for dim=1:size(obj.u,2)     
                    l=obj.degrees(bsis,dim);
                    I_sum=0;
                    for L = 0:l
                        wig_term = Wigner3j([l, l, 2*L], [0, 0, 0]);   
                        I_sum = I_sum + ((4*L + 1) *(wig_term^4));
                    end
                    term2=term2*((2*l+1)^2)*I_sum;
                end
                
                obj.threshs(bsis)=sqrt(term1*term2)/N;
            end
        end

        function C=calcCap(obj,readouts,sample_idxs,basis_idxs,use_bias)
            if sample_idxs==0
                sample_idxs=1:obj.sample_size;
            end
            if basis_idxs==0
                basis_idxs=1:obj.basis_size;
            end

            if use_bias
                X=cat(2,readouts(sample_idxs,:),ones(numel(sample_idxs),1));
            else
                X=readouts(sample_idxs,:);
            end

            N = size(X,1);

            [Q,~]=qr(X,0);
            Xt=Q;

            yy=obj.y(sample_idxs,basis_idxs);
            R=(Xt'*yy);

            Rprod=(1/N)*sum(R.^2,1);
            % ymean=mean(yy.^2,1);
            % C=Rprod./ymean;
            C=Rprod;
        end
     
        function [C_hat,dC_hat]=estCap(obj,X,alg)
            
            C_hat=obj.fitCap(X,1:size(X,1),alg);

            C_hat_1=obj.fitCap(X,1:floor(size(X,1)/2),alg);
            C_hat_2=obj.fitCap(X,floor(size(X,1)/2)+1:size(X,1),alg);

            dC_hat=(1/2)*abs(sum(C_hat_1,"all")-sum(C_hat_2,"all"));

        end
        function C=fitCap(obj,X,sample_idxs,alg)
            if sample_idxs==0
                sample_idxs=1:size(X,1);
            end

            half_len=floor(numel(sample_idxs)/2);
            idxs1=sample_idxs(1:half_len);
            idxs2=sample_idxs(half_len+1:numel(sample_idxs));  

            C_N=obj.calcCap(X,sample_idxs,0,1);       % capacity calculated from whole 
            C_1=obj.calcCap(X,idxs1,0,1);   % capacity calculated with first half
            C_2=obj.calcCap(X,idxs2,0,1);   % capacity calculated with second half

            C_N_avg=mean([C_1;C_2]);   % average half(N/2) estimate

            C=2*C_N-C_N_avg;   %linear fitting

            C(1)=1;     % Capacity of first basis P0=1 due to column of ones

            if alg==1   % algorithm 1: minimum negative as threshold
                zero_idx=C<-min(C);
            elseif alg==2   % algorithm 2: theoretical threshold
                obj.initThresholds(X(sample_idxs,:));
                zero_idx=C_N<obj.threshs;
            end

            C(zero_idx)=0;
         
        end

        function [Cm_arr,Cs_arr]=scanCap(obj,X)
            X=cat(2,X,ones(size(X,1),1));

            max_div=floor(obj.sample_size./(size(X,2)))-1;
            obj.samps_arr=floor(obj.sample_size*(1./(max_div:-1:1)))';
            obj.samps_arr=unique(obj.samps_arr);

            Cm_arr=zeros(length(obj.samps_arr),obj.basis_size); % capacity means
            Cs_arr=zeros(length(obj.samps_arr),obj.basis_size); % capacity stds
            
            for s=1:length(obj.samps_arr)
                if obj.progress_display
                    IPC.dispPerc(s,length(obj.samps_arr));
                end

                N=obj.samps_arr(s);     % number of samples 
                K=floor(obj.sample_size/N); %number of independant subsets
                % K=floor(obj.sample_size/obj.samps_arr(end));

                C_arr=zeros(K,obj.basis_size);

                for k=1:K
                    batch_idx=(k-1)*N+1:k*N;

                    Ct=obj.calcCap(X,batch_idx,0,1);
                    C_arr(k,:)=Ct';
                end
                
                Cm_arr(s,:)=mean(C_arr,1);

                if K>1
                    Cs_arr(s,:)=std(C_arr,1);
                end
            end

        end

        function Cm=Cap_mat(obj,C,K,filename)
            %=====Calculating=====================
            C_tot=sum(C,"all");
            Cm=zeros(obj.max_deg+1);
            % Cm(1,1)=C(1);
            % Cm(1,1)=1;
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
            str = sprintf('$\\frac{%g}{%g}$', round(C_tot,2), K);
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

        function Cd=Cap_deg(obj,C,K,cap_lim,filename)
            % cap_lim=25;
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
            ylim([0,cap_lim]);

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
            str = sprintf('$\\frac{%g}{%g}$', round(C_tot,2), K);
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

        function idx=findBasisIdx(obj,basis_string)
            idx=find(strcmpi(obj.basis_terms,basis_string+" "));
        end

        function orthoPoly3(obj)
            % legendre basis using recurrence relation calculated at data points x
            obj.u_basis=zeros(obj.sample_size,obj.dimn,obj.max_deg+1);
            obj.u_basis(:,:,1) = sqrt(1) * ones(size(obj.u)); 
            obj.u_basis(:,:,2) = sqrt(3) * obj.u;

            for n = 2:obj.max_deg
                IPC.dispPerc(n,obj.max_deg);
                a = sqrt((2*n+1) * (2*n-1)) / n;
                b = ((n-1) / n) * sqrt((2*n+1) / (2*n-3)); 
                obj.u_basis(:,:,n+1) =  a .* obj.u .* obj.u_basis(:,:,n) - b .* obj.u_basis(:,:,n-1);
            end

        end

    end

    
    methods(Static)
        function dispPerc(i,len)
            if(floor(mod(i,len/100))==0)
                percent=i*100/len;
                disp("=== "+percent+" % ===");
            end
        end

        function Z = stars_and_bars(n, k)
            % Generates all ways to distribute n stars and k bars, k bars=k+1 sections

            % Total slots = n stars + k bars
            n_slots= n + k;

            % Get all ways to place k bars among the total slots
            perms = nchoosek(1:n_slots, k);
            num_combos = size(perms, 1);

            Z = zeros(num_combos, k);

            for i = 1:num_combos
                % add start (0) and end (n_slots+1) as boundaries
                bounds= [0, perms(i, :), n_slots + 1];
                for j = 1:k
                    % Stars in section j = gap between consecutive boundaries, minus 1 for the bar itself
                    Z(i, j) = bounds(j + 1) - bounds(j) - 1;
                end
            end
        end

    end

end 





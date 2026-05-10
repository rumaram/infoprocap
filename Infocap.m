classdef Infocap<handle
    properties
        basis_range=[-1,1];     % limits of inner product integral of the hilbert space
        basis_type="poly";
        progress_display=true;

        u;      %inputs.
        y;      %calculated basis functions of inputs u.
        max_deg;
        degrees;
        basis_size;
        sample_size;
        ortho_funcs;
        ortho_u;    % u calculated at all ortho functions
        ortho_syms;
        basis_terms;      
        threshs;
        samps_arr;
        ortho_mat;
    end

    methods
        function obj=Infocap(max_deg)
            obj.max_deg=max_deg;
            if obj.progress_display
                disp("Generating orthogonal basis terms");
            end
            [obj.ortho_funcs,obj.ortho_syms]=Infocap.orthoPoly(obj.max_deg,obj.basis_range,obj.basis_type);
            % obj.orthoPoly2(obj.max_deg);
        end

        function initBasis(obj,u)
            obj.basis_size=nchoosek(obj.max_deg+size(u,2),obj.max_deg);
            
            
            % getting basis
            obj.u=u;
            max_neu=size(u,2);

            obj.degrees=zeros(obj.basis_size,max_neu); 
        
            term_count=0;
            obj.basis_terms={"1"};

            if obj.progress_display
                disp("Calculating basis");
            end
            
            obj.degrees(1,:)=zeros(1,size(obj.degrees,2));
            for d_tot=1:obj.max_deg   % Total degree of the whole polynomial
                
                if(d_tot<max_neu)
                    max_vars=d_tot;
                else
                    max_vars=max_neu;
                end
        
                for vars=1:max_vars    % Number of time steps to use
        
                    d_perms=Infocap.perms_rep(1:d_tot,vars);
                    d_array=[];  % array containing list of individual degress of polynomials
                    for i=1:size(d_perms,1)
                        if(sum(d_perms(i,:))==d_tot)
                            d_array(end+1,:)=d_perms(i,:);
                        end
                    end
        
                    for d=d_array'
                        for win=vars:max_neu    % window= max delay - min dealy +1
        
                            if(vars==1)
                                pos=win;
                            elseif(vars==2)
                                pos=[1,win];
                            else
                                pos=nchoosek(2:win-1,vars-2); % delay of each individual variable
                                pos=[ones(size(pos,1),1),pos,win*ones(size(pos,1),1)];
                                pos=sort(pos,2);
                            end
                            
                                               
                            for p=pos'  % list of delay of each individual variable
                                
                                for neu=0:max_neu   % smallest/starting delay of window
                                        
                                        if(neu+p(end)>max_neu)
                                            break;
                                        end
        
                                        term_count=term_count+1;
        
                                        term="";
                                        for i=1:size(pos,2)
                                            neuron=neu+p(i);
                                            term=term+"P_"+d(i)+"("+"u_"+neuron+")";
                                            if neuron~=0
                                                obj.degrees(term_count+1,neuron)=d(i);
                                            end
                                            
                                        end
                                        
                                        obj.basis_terms{end+1}=term;
                                        if(vars==1)
                                            break;
                                        end
        
                                end
        
                            end
                        end
        
                    end
        
                end
        
            end
            
            obj.basis_size=size(obj.degrees,1);
            obj.basis_terms=string(obj.basis_terms);
            
            % calculating basis for the inputs
            obj.sample_size=size(obj.u,1);
            feature_size=size(obj.u,2);
            
            obj.y=zeros(obj.sample_size,obj.basis_size);

            for basis=1:obj.basis_size
                if obj.progress_display
                    % dispPerc(basis,obj.basis_size);
                end
                prod=ones(obj.sample_size,1);
                
                for feat=1:feature_size
                    if obj.degrees(basis,feat)~=0
                        
                        val=obj.u(:,feat);
                        deg=obj.degrees(basis,feat);
                        ortho_fn=obj.ortho_funcs{deg};

                        lP=ortho_fn(val);   % for ortho fnctions               
                        prod=prod.*lP;
                    end

                end
                obj.y(:,basis)=prod;
            end
            if obj.progress_display
                disp("Basis Initiated");
            end
            
            obj.ortho_mat=(obj.y'*obj.y)./obj.sample_size;

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

            % Xt with svd
            % [U, S, V] = svd(X, 'econ');
            % Xt=U*V';

            % Xt with qr
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

            C_N2=mean([C_1;C_2]);   % half(N/2) sample estimate

            C=2*C_N-C_N2;   %linear fitting

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
                    dispPerc(s,length(obj.samps_arr));
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

        function orthoPoly2(obj, n_max)
            % legendre func handles using recurrence relation
            obj.ortho_funcs=cell(n_max+1,1);
            obj.ortho_funcs{1} = @(x) sqrt(1) * ones(size(x));   % sqrt(2*0+1) = 1
            obj.ortho_funcs{2} = @(x) sqrt(3) * x;

            for n = 2:n_max
                dispPerc(n,n_max);
                a = sqrt((2*n+1) * (2*n-1)) / n;
                b = ((n-1) / n) * sqrt((2*n+1) / (2*n-3));
                f_n   = obj.ortho_funcs{n};      
                f_nm1 = obj.ortho_funcs{n-1};    
                obj.ortho_funcs{n+1} = @(x) a .* x .* f_n(x) - b .* f_nm1(x);
            end
            obj.ortho_funcs(1) = [];
        end

        function orthoPoly3(obj, n_max,x)
            % legendre func handles using recurrence relation
            obj.ortho_u=cell(n_max+1,1);
            obj.ortho_u{1} = sqrt(1) * ones(size(x));   % sqrt(2*0+1) = 1
            obj.ortho_u{2} = sqrt(3) * x;

            for n = 2:n_max
                dispPerc(n,n_max);
                a = sqrt((2*n+1) * (2*n-1)) / n;
                b = ((n-1) / n) * sqrt((2*n+1) / (2*n-3)); 
                obj.ortho_u{n+1} =  a .* x .* obj.ortho_u{n} - b .* obj.ortho_u{n-1};
            end
            obj.ortho_u(1) = [];
        end

    end

    
    methods(Static)

        function C2=Apply_threshold(C,N,T,p)
            % N=number of readouts neurons
            % T=size(obj.u,1);    % T=number of data points
            
            thr=2*chi2inv(1-p,N)/T;
            C2=C;
            for i=1:length(C2)
                if(C2(i)<thr)
                    C2(i)=0;
                end
            end
        end
        
        function y=perms_rep(x,k)
            C = cell(k, 1);             %// Preallocate a cell array
            [C{:}] = ndgrid(x);         %// Create K grids of values
            y = cellfun(@(x){x(:)}, C); %// Convert grids to column vectors
            y = [y{:}];    
        
        end

        function [func_handles,ortho_syms] = orthoPoly(n,range,type)
            a=range(1);
            b=range(2);
            syms t
            assume(t, 'real');
            
            % Initialize the basis with monomials
            basis = sym(zeros(1, n+1));

            if type=="poly"
                % algebraic monomials
                for k = 0:n
                    basis(k+1) = t^k;
                end
            elseif type=="trig"
                % trignometric monomials
                basis(1)=1;
                for k = 1:n
                    basis(2*k)   = cos(k*t);  % Cosine term
                    basis(2*k+1) = sin(k*t);  % Sine term
                end
            end

            % Gram-Schmidt process
            ortho = sym(zeros(1, n+1));
            for k = 0:n
                v = basis(k+1);
                for j = 0:k-1
                    proj = int(ortho(j+1)*v, t, a, b) * ortho(j+1);
                    v = v - proj;
                end
                norm_v = sqrt(int(v^2, t, a, b));
                ortho(k+1) = simplify(v / norm_v);
            end
            ortho=sqrt(2).*ortho;    %normalisation
            ortho_syms=ortho(2:end);
        
            % Step 3: Convert symbolic polynomials to function handles
            % (excluding the first degree function=1)
            func_handles = cell(1, n);
            for k = 1:n
                func_handles{k} = matlabFunction(ortho(k+1), 'Vars', t);
            end
            
        end



    end

end 





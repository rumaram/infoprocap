classdef IPC<handle
    properties
        disp_prog=true; % Display progress in command window
        u;              % inputs.
        u_basis;        % 3d array of first max_deg no. of individual legendre basis functions (not products) calculated at all inputs u
        y;              % calculated product basis functions of inputs u.        
        max_deg;        % highest total degree used to generate basis
        degrees;        % individual degree composition of the product basis terms
        basis_size;     % no of basis terms
        sample_size;    % no of input samples
        dimn;           % input dimensions
        threshs;        % false positive capacity thresholds
        K;              % number of readouts
        basis_names;    % basis term names. Eg: P_1(u_3)
    end

    methods
        function obj=IPC(u,max_deg)
            obj.max_deg=max_deg;
            obj.u=u;
            obj.sample_size=size(u,1);
            obj.dimn=size(u,2);
         
            obj.degrees=infoprocap.Utils.stars_and_bars(obj.max_deg,obj.dimn); 
            % distribution of degrees of product basis is equivalent to a stars and bars Combinatorial problem
            
            obj.basis_size=size(obj.degrees,1);
            obj.y=ones(obj.sample_size,obj.basis_size);

            obj.threshs=zeros(1,obj.basis_size);

            %==Calculating legendre basis using recursive relations calculated at all input samples==
            % Normalization: 0.5*integ_{-1}^{1} Pm(x)*Pn(x) dx=delta_{mn}
            obj.u_basis=zeros(obj.sample_size,obj.dimn,obj.max_deg+1);  % samples x dimensions x degrees
            obj.u_basis(:,:,1) = sqrt(1) * ones(size(obj.u)); 
            obj.u_basis(:,:,2) = sqrt(3) * obj.u;

            for n = 2:obj.max_deg            
                a = sqrt((2*n+1) * (2*n-1)) / n;
                b = ((n-1) / n) * sqrt((2*n+1) / (2*n-3)); 
                obj.u_basis(:,:,n+1) =  a .* obj.u .* obj.u_basis(:,:,n) - b .* obj.u_basis(:,:,n-1);
            end

            %==Calculating product basis terms and basis names==
            obj.basis_names=strings(1,obj.basis_size);
            if obj.disp_prog
                disp("Basis initialization Progress: ");
                disp('     ');
            end
            for basis=1:obj.basis_size
                if obj.disp_prog
                    infoprocap.Utils.dispPerc(basis,obj.basis_size);
                end
                obj.basis_names(basis)="";
                for q=1:obj.dimn             
                    deg=obj.degrees(basis,q);
                    obj.y(:,basis)=obj.y(:,basis).*obj.u_basis(:,q,deg+1);
                    if deg~=0
                        obj.basis_names(basis)=obj.basis_names(basis)+"P_"+deg+"(u_"+q+")";
                    end
                    obj.basis_names(1)="1";
                end             
            end
            
            if obj.disp_prog
                disp(' ');
                disp("Basis Initialized");
            end

        end

        function initThresholds(obj,X)
            % Initialize thresholds for Algorithm 2

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
                        wig_term = infoprocap.Utils.W3j(l, 2*L);   
                        I_sum = I_sum + ((4*L + 1) *(wig_term^4));
                    end
                    term2=term2*((2*l+1)^2)*I_sum;
                end
                
                obj.threshs(bsis)=sqrt(term1*term2)/N;
            end
        end

        function C=calcCap(obj,readouts,sample_idxs,basis_idxs,use_bias)
            % Calculates raw capacities
            % use_bias: Add an extra column of ones to readouts

            arguments
                obj 
                readouts 
                sample_idxs =1:obj.sample_size
                basis_idxs =1:obj.basis_size
                use_bias =1
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
            % C=Rprod./ymean;   % normalization ignored since E[y^2]=1
            C=Rprod;
        end

        function C=estCap(obj,X,alg,sample_idxs)
            % Fit capacities based on the selected algorithm
            % Algorithm 1: Theoretical threshold, 2: Minimum negative threshold

            arguments
                obj 
                X 
                alg =1
                sample_idxs=1:size(X,1);         
            end
            obj.K=size(X,2);

            half_len=floor(numel(sample_idxs)/2);
            idxs1=sample_idxs(1:half_len);
            idxs2=sample_idxs(half_len+1:numel(sample_idxs));  

            C_N=obj.calcCap(X,sample_idxs);       % capacity calculated from whole 
            C_1=obj.calcCap(X,idxs1);   % capacity calculated with first half
            C_2=obj.calcCap(X,idxs2);   % capacity calculated with second half

            C_N_avg=mean([C_1;C_2]);   % average half(N/2) estimate

            C=2*C_N-C_N_avg;   %linear fitting

            C(1)=1;     % Capacity of first basis P0=1 due to column of ones

            if alg==1   
                obj.initThresholds(X(sample_idxs,:));
                zero_idx=C_N<obj.threshs;   % threshold compared with raw capacities.   
            elseif alg==2  
                zero_idx=C<-min(C);         % threshold compared with fitted capacities.
            end

            C(zero_idx)=0;
            C(C<0)=0;
         
        end

        function [samps_arr,Cm_arr]=scanCap(obj,X)
            % scans capacities changing number of samples. 
            % Useful to observe asymptotic form of capacities.

            X=cat(2,X,ones(size(X,1),1));

            max_div=floor(obj.sample_size./(size(X,2)))-1;
            samps_arr=floor(obj.sample_size*(1./(max_div:-1:1)))';
            samps_arr=unique(samps_arr);

            Cm_arr=zeros(length(samps_arr),obj.basis_size); % capacity means
            
            if obj.disp_prog
                disp("Scan Progress: ");
                disp('     ');
            end
            for s=1:length(samps_arr)
                     
                if obj.disp_prog
                    infoprocap.Utils.dispPerc(s,length(samps_arr));
                end
                

                n_samp=samps_arr(s);     % number of samples 
                n_parts=floor(obj.sample_size/n_samp); %number of independent partitions
                C_arr=zeros(n_parts,obj.basis_size);

                for p=1:n_parts
                    batch_idx=(p-1)*n_samp+1:p*n_samp;

                    Ct=obj.calcCap(X,batch_idx);
                    C_arr(p,:)=Ct';
                end
                
                Cm_arr(s,:)=mean(C_arr,1);
            end
            disp(' ');
            disp("Capacity scan finished");
        end

        function idx=nameToidx(obj,name)    % converts basis name to corresponding index
            idx=find(obj.basis_names,name);
        end

        function name=idxToname(obj,idx)    % converts basis index to name of the basis term
            name=obj.basis_names(idx);
        end

        function exportCSV(obj, C, filename)
        % Exports IPC results to a CSV file.
        % Thresholds are all zeros if initThresholds has not been called, for eg. if Algorithm 2 is used
       
            fid = fopen(filename, 'w');
            if fid == -1
                error('IPC:exportCSV:cannotOpenFile', ...
                      'Could not open file for writing: %s', filename);
            end
        
            % --- Header block ---
             fprintf(fid, '---IPC results--- \n');
            fprintf(fid, '\n');
            fprintf(fid, 'Total Capacity,%.2f\n', sum(C, "all"));
            fprintf(fid, 'Number of nonzero capacities,%d\n',  nnz(C));
            fprintf(fid, '\n');
            fprintf(fid, 'Number of input samples,%d\n',  obj.sample_size);
            fprintf(fid, 'Number of readouts,%d\n',        obj.K);
            fprintf(fid, 'Number of input dimensions,%d\n', obj.dimn);
            fprintf(fid, 'Maximum total degree,%d\n',       obj.max_deg);
            fprintf(fid, 'Number of basis functions,%d\n',  obj.basis_size);      
            fprintf(fid, '\n');
        
            % --- Column titles ---
            fprintf(fid, 'Basis Index,Basis Name,Capacity,Threshold\n');
        
            % --- One row per basis function ---
            for i = 1:obj.basis_size      
                fprintf(fid, '%d,%s,%.6f,%.6f\n', i,obj.basis_names(i), C(i),obj.threshs(i));
            end
        
            fclose(fid);
            fprintf('IPC results saved to: %s\n', filename);
        end

    end

end 





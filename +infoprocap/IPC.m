classdef IPC<handle
    properties
        progress_display=true;
        u;              % inputs.
        u_basis;        % cell array of first max_deg no. of individual legendre basis functions (not products) calculated at input u
        y;              % calculated product basis functions of inputs u.        
        max_deg;        % highest total degree used to generate basis
        degrees;        % individual degree composition of the product basis
        basis_size;     % no of basis terms
        sample_size;    % no of input samples
        dimn;           % input dimension
        threshs;        % false positive capacity thresholds
        K;              % number of readouts
    end

    methods
        function obj=IPC(u,max_deg)
            obj.max_deg=max_deg;
            obj.u=u;
            obj.sample_size=size(u,1);
            obj.dimn=size(u,2);

            % distribution of basis is equivalent to a stars and bars Combinatorial problem
            obj.degrees=infoprocap.Utils.stars_and_bars(obj.max_deg,obj.dimn);  
            
            obj.basis_size=size(obj.degrees,1);
            obj.y=ones(obj.sample_size,obj.basis_size);

            %==Calculating legendre basis using recursive relations==
            % Normalization: 0.5*integ_{-1}^{1} Pm(x)*Pn(x) dx=delta_{mn}
            obj.u_basis=zeros(obj.sample_size,obj.dimn,obj.max_deg+1);  % samples x dimensions x degrees
            obj.u_basis(:,:,1) = sqrt(1) * ones(size(obj.u)); 
            obj.u_basis(:,:,2) = sqrt(3) * obj.u;

            for n = 2:obj.max_deg            
                a = sqrt((2*n+1) * (2*n-1)) / n;
                b = ((n-1) / n) * sqrt((2*n+1) / (2*n-3)); 
                obj.u_basis(:,:,n+1) =  a .* obj.u .* obj.u_basis(:,:,n) - b .* obj.u_basis(:,:,n-1);
            end

            %==Calculating product basis terms==
            if obj.progress_display
                disp("Basis initialization Progress: ");
                disp('     ');
            end
            for basis=1:obj.basis_size
                if obj.progress_display
                    infoprocap.Utils.dispPerc(basis,obj.basis_size);
                end
                for q=1:obj.dimn             
                    deg=obj.degrees(basis,q);
                    obj.y(:,basis)=obj.y(:,basis).*obj.u_basis(:,q,deg+1);
                end             
            end
            
            if obj.progress_display
                disp(' ');
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
            % C=Rprod./ymean;
            C=Rprod;
        end
     
        function [C_hat,dC_hat]=estCap(obj,X,alg)
            % C_hat = estimated capacity based on the selected algorithm
            % dC_hat = an approximate uncertainty estimate for capacities by fitting capacities on two halves

            obj.K=size(X,2);
            
            C_hat=obj.fitCap(X,1:size(X,1),alg);

            C_hat_1=obj.fitCap(X,1:floor(size(X,1)/2),alg);
            C_hat_2=obj.fitCap(X,floor(size(X,1)/2)+1:size(X,1),alg);

            dC_hat=(1/2)*abs(sum(C_hat_1,"all")-sum(C_hat_2,"all"));

        end

        function C=fitCap(obj,X,sample_idxs,alg)
            % fit capacities based on the selected algorithm
            arguments
                obj 
                X 
                sample_idxs=1:size(X,1);
                alg =1
            end

            half_len=floor(numel(sample_idxs)/2);
            idxs1=sample_idxs(1:half_len);
            idxs2=sample_idxs(half_len+1:numel(sample_idxs));  

            C_N=obj.calcCap(X,sample_idxs);       % capacity calculated from whole 
            C_1=obj.calcCap(X,idxs1);   % capacity calculated with first half
            C_2=obj.calcCap(X,idxs2);   % capacity calculated with second half

            C_N_avg=mean([C_1;C_2]);   % average half(N/2) estimate

            C=2*C_N-C_N_avg;   %linear fitting

            C(1)=1;     % Capacity of first basis P0=1 due to column of ones

            if alg==1   % algorithm 1: theoretical threshold
                obj.initThresholds(X(sample_idxs,:));
                zero_idx=C_N<obj.threshs;          
            elseif alg==2   % algorithm 2: minimum negative as threshold
                zero_idx=C<-min(C);
            end

            C(zero_idx)=0;
         
        end

        function [samps_arr,Cm_arr]=scanCap(obj,X)
            % scans capacities changing number of samples. Useful to observe asymptotic form
            X=cat(2,X,ones(size(X,1),1));

            max_div=floor(obj.sample_size./(size(X,2)))-1;
            samps_arr=floor(obj.sample_size*(1./(max_div:-1:1)))';
            samps_arr=unique(samps_arr);

            Cm_arr=zeros(length(samps_arr),obj.basis_size); % capacity means
            
            for s=1:length(samps_arr)
                
                disp("Scan Progress");
                disp('     ');
                if obj.progress_display
                    infoprocap.Utils.dispPerc(s,length(samps_arr));
                end
                disp(' ');

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

        end

    end

end 





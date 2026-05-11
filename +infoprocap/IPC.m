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

            obj.degrees=infoprocap.Utils.stars_and_bars(obj.max_deg,obj.dimn);  % dsitribution of basis is equivalent to a stars and bars problem
            
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
                    infoprocap.Utils.dispPerc(s,length(obj.samps_arr));
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

        function idx=findBasisIdx(obj,basis_string)
            idx=find(strcmpi(obj.basis_terms,basis_string+" "));
        end

        function orthoPoly3(obj)
            % legendre basis using recurrence relation calculated at data points x
            obj.u_basis=zeros(obj.sample_size,obj.dimn,obj.max_deg+1);
            obj.u_basis(:,:,1) = sqrt(1) * ones(size(obj.u)); 
            obj.u_basis(:,:,2) = sqrt(3) * obj.u;

            for n = 2:obj.max_deg
                infoprocap.Utils.dispPerc(n,obj.max_deg);
                a = sqrt((2*n+1) * (2*n-1)) / n;
                b = ((n-1) / n) * sqrt((2*n+1) / (2*n-3)); 
                obj.u_basis(:,:,n+1) =  a .* obj.u .* obj.u_basis(:,:,n) - b .* obj.u_basis(:,:,n-1);
            end

        end

    end

end 





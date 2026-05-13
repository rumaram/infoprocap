classdef Phot_sys<handle

    properties

        N           = 2^14;
        dL     = 50e-2;    % Small segment for each step of split-step fourier method
        L=5;    % Total propagation length
        T=1000e-12;
        dis_save=50e-2;  % distance in between to save results 
        P_avg = 6.6;      % Avg_power in dBm at the output

        wav_centr=1550e-9;
        wav_span_encode=2.5e-9;
        wav_span_meas=3.5e-9;  % wavelength span to take measurements
        rep_rate=10e6;
        
        waveshaper_res_hz=12e9;    % frequency resolution of waveshaper in Hz
        osa_res_nm=0.05;
       
        bandW               = 0.6e-9;      % bandwidth in wavelength   
        gamma       = 1.2e-3;
        beta2=-2.3e-26;

        % Declarations
        t
        L_d;
        L_nl;
        pulseW;     %fwhm
        peakP; 
        span_encode;
        span_meas;
        E;
        Ef; 
        expD; % Dispersion operator
        wavs;    % to store wavelengths array centered at wav_centr
        wavs_meas;  % wavelengths of measured span

        readouts;   % Stores spectral intensity readouts at each length step
        X;
        w;
        dw;
    end



    properties(Constant=true)
        c=2.99792458e8;
    end

    
    methods
        function obj=Phot_sys()     
            obj.updateParams();
        end

        function updateParams(obj)           
            rng(3);

            dt=obj.T/obj.N;
            obj.t=(-obj.N/2:obj.N/2-1)*dt;
            obj.dw=2*pi/obj.T;
            obj.w=(-obj.N/2:obj.N/2-1)*obj.dw;

            obj.pulseW=obj.bandW2pulseW(obj.bandW,"sech2"); 
            obj.peakP=obj.avg2peak(obj.P_avg,"sech2");
            t0=obj.pulseW/1.763;

            obj.E=sqrt(obj.peakP)* sech(obj.t/t0);

            freq0=obj.c/obj.wav_centr;
            freqs=(obj.w./(2*pi))+freq0;
            obj.wavs=obj.c./freqs;
            obj.wavs=obj.wavs(end:-1:1);

            [obj.span_encode,~]=obj.getSpanidx(obj.wav_span_encode);
            [obj.span_meas,obj.wavs_meas]=obj.getSpanidx(obj.wav_span_meas);          

            obj.Ef=fftshift(fft(obj.E));

            obj.L_d=(obj.pulseW/1.76)^2/abs(obj.beta2);
            obj.L_nl=1/(obj.gamma*obj.peakP);
            
            D=(1i*obj.beta2*(obj.w.^2)/2);
            obj.expD=exp(D*obj.dL/2);
        end
        
        function run(obj,X)

            obj.X=X;
            sample_size=size(X,1);

            len_save=floor(obj.L/obj.dis_save);    % number of lengths to save

            R=zeros(len_save+1,length(obj.span_meas),sample_size);
            fourrier_norm=(obj.N).^2;
            
            disp("Photonic system run Progress:");
            disp('     ');
            for i=1:sample_size           
                infoprocap.Utils.dispPerc(i,sample_size);
                modulated_field=obj.waveShape(X(i,:));
                propagated_field=obj.propLight(modulated_field,obj.L);
                propagated_spectrum= fftshift(fft(propagated_field,[],2),2);
                output=abs(propagated_spectrum(:,obj.span_meas)).^2./(fourrier_norm);

                R(:,:,i)=output;  
            end
            disp(' ');

            obj.readouts=R;
            disp("Photonic system run finished");
 
        end


        function A_prop=propLight(obj,A,distance)  %Split-step fourrier propagation
            len_prop=length(0:obj.dL:distance);    % number of lengths to propagate
            len_save=floor(distance/obj.dis_save);    % number of lengths to save
            dis_bin=floor(obj.dis_save/obj.dL);     % bin size corresponding to save distance
            
            A_prop=zeros(len_save+1,length(A));
            Ah=A;
            A_prop(1,:)=Ah;

            dis=2;
            for i=2:len_prop
                Ah=Dispers(Nonlin(Dispers(Ah,obj.expD)),obj.expD);

                if(mod(i,dis_bin)==0)
                    A_prop(dis,:)=Ah;
                    dis=dis+1;
                end
            end

            function Ah=Dispers(A,expD2) % exponential dispersion function, input & output in time domain
                product_centered=expD2.*fftshift(fft(A));% expD is already centered. fftshift make the 2nd term centered to match.
                Ah=ifft(ifftshift(product_centered));%ifftshift uncenter the product
            end
            
            function Ah=Nonlin(A) % exponential nonlinearity function, input & output in time domain
                I=abs(A).^2;    %Intensity
                NL=1i*obj.gamma*I;                
                Ah=exp(NL*obj.dL).*A;
            end
        
        end

        function E_modulated=waveShape(obj,x)
            N_span_encode=length(obj.span_encode);
            fill_idx=obj.fillBin(N_span_encode,length(x));
            
            mask=ones(1,N_span_encode);
            
            for i=1:N_span_encode
                mask(i)=x(fill_idx(i));
            end
            mask=sign(mask).*sqrt(abs(mask));

            % making mask same size as E by padding ones
            mask2=ones(size(obj.E));    
            mask2(obj.span_encode)=mask;
            mask=mask2;
            %===
            mask_amp=abs(mask);

            mask_ph=sign(mask);
            mask_ph(mask_ph>=0)=0;
            mask_ph(mask_ph<0)=pi;

            phase_noise=normrnd(0,0.15e-2*2*pi,size(mask_ph));
            mask_ph=mask_ph+phase_noise;

            mask=mask_amp.*exp(1i*mask_ph);

            span_padded=obj.getSpanidx(obj.wav_span_encode+1e-9); % higher pad span to reduce edge effects
            span_padded2=obj.getSpanidx(obj.wav_span_encode+0.5e-9);  %lower pad span to include overflow outside encoding span
            
            mask2=mask;
            mask2(span_padded)=obj.applyIF(mask(span_padded),obj.waveshaper_res_hz,"flattop");
            mask(span_padded2)=mask2(span_padded2);

            Ef2=fftshift(fft(obj.E));
            Ef2=Ef2.*mask;
            E_modulated=ifft(ifftshift(Ef2));      
            
        end

        % ====Utility functions=============%

        function R=prepReadouts(obj,dis)  % Prepare readouts at a particular fiber length
            dis_idx=floor(dis/obj.dis_save)+1;
            readouts2=permute(obj.readouts,[3,2,1]);
            R=readouts2(:,:,dis_idx);
            R=R(:,end:-1:1);

            wavs_interp=obj.wavs_meas(1):obj.osa_res_nm*1e-9:obj.wavs_meas(end);    % matching OSA resolution
            R=interp1(obj.wavs_meas, R', wavs_interp, 'linear');   
            R=R';
        end

        function dt=bandW2pulseW(obj,dwav,shape)   % converts pulse width in wavelength to time
            % spec_wid in wavelength
            if shape=="sech2"
                tbp=0.315;
            end
            if shape=="gaussian"
                tbp=0.441;
            end

            f0=obj.c/obj.wav_centr;
            wav1=obj.wav_centr+dwav;
            f1=obj.c/wav1;
            df=f0-f1;
            
            dt=tbp/df;
        end

        function P_peak=avg2peak(obj,P_avg,mode)
            % P_avg in dbm, P_peak in W, fwhm in seconds
            P_avg_w=10^(P_avg/10)/1000; % P_avg in watt
            if (mode=="sech2")
                P_peak=0.8815*P_avg_w/(obj.rep_rate*obj.pulseW);
            elseif(mode=="gaussian")
                P_peak=0.94*P_avg_w/(obj.rep_rate*obj.pulseW);
            end
        end

        function C=fillBin(obj,lenA,lenB)   % return indices of B that fills A by uniformly repeating B. (lenA>lenB)
    
            C=zeros(1,lenA);
        
            bin_width_low=floor(lenA/lenB);
            bin_width_high=ceil(lenA/lenB);
            
            frac=(lenA/lenB)-floor(lenA/lenB);  %amount left behind when you use bin_width_low.
            decimal_part=0; %stores accumulated decimal parts left behind due to rounding of bin_width
            idx_start=1; 
            idx_B=1;
        
            while idx_B<=lenB
                if decimal_part>1   %when decimal part become>1 u use bin_width_high to compensate for lost parts
                    bin_width=bin_width_high;
                    decimal_part=decimal_part-1+frac;
                else
                    bin_width=bin_width_low;
                    decimal_part=decimal_part+frac;
                end
                C(idx_start:idx_start+bin_width-1)=idx_B;
                idx_start=idx_start+bin_width;
                idx_B=idx_B+1;
            end
            C(C==0)=idx_B-1;
    
        end

        function [span,wavs]=getSpanidx(obj,wav_span)
            wav_left=obj.wav_centr-wav_span./2;
            wav_right=obj.wav_centr+wav_span./2;
            freq_left=obj.c/wav_right;
            freq_right=obj.c/wav_left;

            freq0=obj.c/obj.wav_centr;
            freqs=(obj.w./(2*pi))+freq0;

            [~,left_idx]=min(abs(freqs-freq_left));
            [~,right_idx]=min(abs(freqs-freq_right));

            span=left_idx:right_idx;

            wavs=obj.c./freqs(span);
            wavs=wavs(end:-1:1);
        end

        function phi= getNLphase(obj,dis) % Get nonlinear phase at a length in multiple of pi
            phi=(1/obj.L_nl)*dis/pi;
        end

        function B=applyIF(obj,A,filter_res,type)    % instrument filter (gaussian convolution)
            %filter_res= resolution of the filter in Hz
            
            df=obj.dw/(2*pi);
            N_filter=ceil(filter_res/df);     %number of samples in filter resolution
            N_sig = N_filter / (2*sqrt(2*log(2))); % standard deviation in units of no. of samples 
            
            if N_filter<1
                disp("Incompatible inputs");
                return;
            end
            
            if type=="gaussian"
                % Build Gaussian kernel
                half_width = ceil(4*N_sig);     % cover ±4σ
                x=-half_width:half_width;
                kernel=exp(-0.5*(x/N_sig).^2);
                kernel = kernel / sum(kernel);   % normalize to unit area

            elseif type=="flattop"
                kernel = ones(1, N_filter) / N_filter; 

            elseif type=="tukey"    % flat top with curved edges
                flat_fraction=0.7;  % fraction of window that is flat
                kernel = tukeywin(N_filter, 1-flat_fraction);    % tukey kernal
                kernel = kernel / sum(kernel);
            end
            B=conv(A,kernel,'same');
            end

    end

end


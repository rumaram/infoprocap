classdef Utils
    methods(Static)

        function dispPerc(i,len)
            if(floor(mod(i,len/100))==0)
                fprintf('\b\b\b\b\b\b%05.2f%%', i/len*100);
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

        function W=W3j(j1,j3)
            % Special Wigener3j function when j1=j2, and j1,j2,j3 are non negative integers
            if j3 > 2*j1 || rem(j3,2)
                W=0;
                return;
            end

            t1 = j1-j3;
            t3 = 2*j1-j3;

            t = max( 0, t1 ) : min( t3, j1 );

            fac_t = -( gammaln(t+1) + 2*gammaln(t-t1+1) + gammaln(t3-t+1) + 2*gammaln(j1-t+1) );
            fac_j = 0.5 * ( -gammaln(2*j1+j3+2) + gammaln(2*j1-j3+1)+ 4*gammaln(j3+1) + 4*gammaln(j1+1) );

            W = sum( (-1).^t .* exp(fac_t+fac_j ) );

        end

    end
end
classdef Utils
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
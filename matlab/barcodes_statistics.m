% This program reads multiple files, each containing counts for a specific 
% barcode. The program estimates the SHAPE reactivities for each bar code 
% and compares them to the estimated reactivities for the WT molecule. 
close all
clear all

basedir = '';

% ## SEQUENCE PROPERTIES ##
% The optimization flag options are: 
% 1 = set negatives to zero and normalize, 
% 2 = apply numerical optimization.
optim_flag = 2;  
% This is the number of iterations to be used by the numerical optimization
% procedure.
num_iter = 20;  % this is the max number of iterations allowed.
accuracy = 0.0000001;  % this is the desired optimization accuracy.
n = 172; % number of nucleotides to be read from the data file
analyzed_sequence_start_index = 2;  % 5' (top)
analyzed_sequence_end_index = 172;  % 3' (bottom)

nucleotides = ['A' 'G' 'C' 'T'];
counts_data = zeros(n, 4); 
signal_data = zeros(n, 4); 
WT_Theta = zeros(n-1, 2);
ML_Theta = zeros(n-1, 256);
ML_Gamma =  zeros(n-1, 256);
ML_c = zeros(1, 256);
intermediate_c = zeros(1, 256);
total_counts = zeros(2, 256);
JS_divergences = zeros(1, 256);
p_0_hats = zeros(1, 256);
sum_negatives = zeros(1, 256);
num_negatives = zeros(1, 256);
barcodes = char(256, 4);
counter = 0;

% ## RETRIEVE WT STATISTICS ## 
fid = fopen(strcat(basedir, 'target_WT.adducts'), 'r');
% Skip the contents of the first line since it pertains to complete 
% fragments (we retrieve to advance the file pointer). 
first_line = fgetl(fid);
for i=1:(n-1)
    % there are 2 columns per line: nucleotide ID and Theta value
    WT_Theta(i,:) = fscanf(fid, '%g %g', [1 2]);  
end
status = fclose(fid);

WT_Theta = flipud(WT_Theta);

% ## LOOP OVER ALL BARCODE DATA ##
for letter1=1:4
    for letter2=1:4
        for letter3=1:4
            for letter4=1:4
                counter = counter + 1;
                barcode = strcat(nucleotides(letter1), nucleotides(letter2), nucleotides(letter3), nucleotides(letter4));
                                    
                % ## READ COUNTS ##
                % File format:   Index     Base    (+)-Channel_Count     (-)-Channel_Count
                % Bases are ordered from 5' to 3', top to bottom. 
                fid = fopen(strcat(basedir, 'target_', barcode, '.adducts'), 'r');
                if (fid > 2)
                    barcodes(counter,1:4) = barcode;
                else
                    barcodes(counter,1:4) = 'NNNN';
                end
                % get the first line to advance the file pointer to the data section
                first_line = fgetl(fid);
                for i=1:n
                     % there are 4 columns per line, where the last two contain the counts
                     counts_data(i,:) = fscanf(fid, '%d %s\t %d %d', [1 4]);  
                end
                status = fclose(fid);

                % Organize data by direction of signal decay (i.e., from position 1 to 
                % position n).
                signal_data = flipud(counts_data);
                total_counts(counter, 1) = sum(counts_data(:,3));  % (+) channel
                total_counts(counter, 2) = sum(counts_data(:,4));  % (-) channel

                % ## PRE-PROCESS DATA ##
                % Here, we omit from analysis all nucleotides that had zero counts in both
                % channels (X_k = Y_k = 0). Their reactivities are automatically set to 0. 
                % This way, we make sure each position k in the resulting sequence has at
                % least one non-zero count. This is required for the feasibility and 
                % stability of the numerical optimization procedure. Later, we re-assemble 
                % the entire sequence. 
                analyzed_ind = sort(union(find(signal_data(:, 3) > 0), find(signal_data(:, 4) > 0)));
                analyzed_data = zeros(length(analyzed_ind), 4);
                analyzed_data = signal_data(analyzed_ind, :);
                
                % Locate all positions for which X_k=0 or Y_k=0, and set the counts to 1.
                % This correction is applied for numerical stability, but can be omitted
                % later!
                zero_Ys = find(analyzed_data(:,4) == 0);
                analyzed_data(zero_Ys, 4) = 1;
                zero_Xs = find(analyzed_data(:,3) == 0);
                analyzed_data(zero_Xs, 3) = 1;
                
                % This is the length of the analyzed sequence, including the complete
                % fragments counts, which are at the end.
                N = length(analyzed_ind);
                opt_Theta = zeros(N-1, 1);
                opt_Gamma = zeros(N-1, 1);
                                
                % ## COMPUTE MLE ##
                [opt_Theta, opt_Gamma, ML_c(counter), intermediate_c(counter), sum_negatives(counter), num_negatives(counter), p_0_hats(counter), status_flag] = MLE_one_molecule(analyzed_data, optim_flag, num_iter, accuracy);
                 
                % ## RE-ASSEMBLE ORIGINAL SEQUENCE ##
                ML_Theta(analyzed_ind(1:N-1),counter) = opt_Theta;
                ML_Gamma(analyzed_ind(1:N-1),counter) = opt_Gamma;
                                
                % ## COMPUTE JS DISTANCE TO WT DISTRIBUTION ## 
                if (status_flag > 0)
                    JS_divergences(counter) = entropy(0.5*(ML_Theta(:,counter)+WT_Theta(:,2)))- 0.5*entropy(ML_Theta(:,counter)) - 0.5*entropy(WT_Theta(:,2));
                else  % this indicates that the data set was not valid (e.g., had zero complete fragments or major inconsistencies)
                    JS_divergences(counter) = -999;
                end
            end
        end
    end
end

save workspace



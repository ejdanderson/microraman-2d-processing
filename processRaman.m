addpath(genpath(strcat(fileparts(which(mfilename)), '/subtightplot')));
data = 'F:\Dropbox\Dropbox\Raman\20180305-ascii\halfmicronscan';
%data = '/Users/evan/Dropbox/Raman/20180305-ascii/halfmicronscan';

delimiterIn = '\t';
headerlinesIn = 3;
A = importdata(data, delimiterIn, headerlinesIn);
num_steps_x = str2num(A.textdata{2});
num_steps_y = str2num(A.textdata{1});

matSize = num_steps_x * num_steps_y;

showAllSpectraFig = false;

Spectrum = zeros(matSize, 1340);

% Peak to integrate over (in pixels) and the title
peaks = {
    %35 75 'PO'; ...
    %140 160 'PO'; ...
    400 440 'PO4'};
    %460 485 'PO'; ...
    %490 520 'CO3'; ...
    %615 660 'Amide III'; ...
    %770 810 'CH2 / Amide II'; ...
    %930 990 'Amide I'};

   
grids = zeros(num_steps_x, num_steps_y, length(peaks));

%{
This is the threshold for a single pixel to increase by before being 
subtracted (set to the previous pixel value). This removes cosmic rays to a 
degree; occasionally they will be two pixels wide.
%}
cosmicThreshold = 300;

% Create a matrix with the spectra for each point as rows 
for i = 0:num_steps_y-1
    for j=1:num_steps_x
        k = j+num_steps_x*i;
        l = k;%k*2-1;
        Spectrum(k,:)=A.data(l,:);
    end
end

% Remove bad rows such as focusing on an imperfection and creating 100%
% fluorescence
%Spectrum(207, :) = Spectrum(206, :);
%Spectrum(400, :) = Spectrum(399, :);

rowsSpectrum = size(Spectrum, 1);
colsSpectrum = size(Spectrum, 2);

for j=1:colsSpectrum
    minVal = 999999;
    for i=1:rowsSpectrum
        if i ~= 1 && i ~= rowsSpectrum && j ~= 1 && j ~= colsSpectrum
            lastVal = Spectrum(i - 1, j - 1);
            thisVal = Spectrum(i, j);
            nextVal = Spectrum(i + 1, j + 1);  
            if thisVal - lastVal > cosmicThreshold && thisVal - nextVal > cosmicThreshold
                Spectrum(i, j) = lastVal;
            end
        end
    end
    %Spectrum(i,:) = Spectrum(i,:) / max(Spectrum(i,:));
    %Spectrum(i, :) = (zeros(1, colsSpectrum) + minVal);
end

spectra_out = zeros(num_steps_x, size(Spectrum, 2));
testSpectrum = zeros(size(Spectrum,1), 41);
for peak_num=1:size(peaks, 1)
    for x = 0:num_steps_x-1
        for y = 1:num_steps_y
            index = (num_steps_y * x) + y;
            
            if mod(x, 2) == 0
                y_pos = num_steps_y - y + 1;
            else
                y_pos = y;
            end
            grids(x+1, y_pos, peak_num) = sum(Spectrum(index, peaks{peak_num, 1}:peaks{peak_num, 2}));
            if peak_num == 3
               testSpectrum(index, :) = Spectrum(index, peaks{peak_num, 1}:peaks{peak_num, 2});
               testSpectrum(index, :) = testSpectrum(index, :) - (Spectrum(index, peaks{peak_num, 1})+Spectrum(index, peaks{peak_num, 1}))/4;
            end
            
            %This attempts to normalize the spectra by removing the
            %baseline - i.e. subtract the average of the starting and
            %ending intensity
            grids(x+1, y_pos, peak_num) = grids(x+1, y_pos, peak_num) - (Spectrum(index, peaks{peak_num, 1})+Spectrum(index, peaks{peak_num, 1}))/2 * (peaks{peak_num, 2} - peaks{peak_num, 1});
            %if y_pos == 5
            spectra_out(x+1,:) = spectra_out(x+1,:) + Spectrum(index,:);
           % end
                
            %%sum(Spectrum(index, peaks{peak_num, 1}:peaks{peak_num, 2}))
        end
    end
end

%{
Remove cosmic ray values. This works by checking if adjacent pixels to the
current one change by the `cosmicThreshold` value in a peak like manner.  

e.g. pixels p1, p2, p3 have values of 0, 10000, and 50 respectivly. 
%}


figure
subplot_size = ceil(sqrt(size(peaks,1)));
if size(peaks, 1) > 1
    for peak_num=1:size(peaks, 1)   
        subtightplot(3, 3, peak_num, [0.01, 0.03, 0.03])
        imagesc(transpose(grids(:, :, peak_num)))
        pbaspect([1 num_steps_y/num_steps_x 1])
        colorbar
        title(peaks(peak_num, 3))
    end
else
    imagesc(transpose(grids(:, :, 1)))
    pbaspect([1 num_steps_y/num_steps_x 1])
    colorbar
    caxis([0 1.5e5]);
    title(peaks(peak_num, 3))
end
if showAllSpectraFig
    figure
    zz = transpose(Spectrum);
    numRows = size(Spectrum, 1);
    numCols = size(Spectrum, 2);

    xx = zeros(numRows, numCols);
    yy = zeros(numRows, numCols);


    for i=1:numRows
        for j=1:numCols
            xx(i, j) = j;
            yy(i, j) = i;
        end
    end
    plot3(xx, yy, zz, 'Color', [0.2 0.2 0.8])

    xlabel('pixel')
    ylabel('sample #')
    zlabel('Intensity (arb. units)')
end

dlmwrite('F:\Dropbox\Dropbox\Raman\20180305-ascii\tip-halfmicron-mat', spectra_out);

x_vals = importdata('F:\Dropbox\Dropbox\Teeth\scans\xvals.txt');
figure
plot(1:1340, spectra_out)

if false
    for peak_num=1:size(peaks, 1)
        figure
        plot(grids(:, :, peak_num))
    end
end

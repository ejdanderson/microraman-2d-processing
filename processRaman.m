addpath(genpath(strcat(fileparts(which(mfilename)), '/subtightplot')));
data = '/Users/evan/Dropbox/Teeth/20171017/toothB-scan.txt';

delimiterIn = '\t';
headerlinesIn = 3;
A = importdata(data, delimiterIn, headerlinesIn);
num_steps_x = str2num(A.textdata{1});
num_steps_y = str2num(A.textdata{2});

matSize = num_steps_x * num_steps_y;

showAllSpectraFig = false;

% Peak to integrate over (in pixels) and the title
peaks = {48 70 'PO'; ...
    140 160 'PO'; ...
    400 440 'PO4'; ...
    460 490 'PO'; ...
    500 520 'CO3'; ...
    610 615 'Amide III'; ...
    775 805 'CH2 / Amide II'; ...
    883 900 '?'; ...
    930 990 'Amide I'};

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
        l = (i*((num_steps_x+1)*2-1)+j*2)-1;
        Spectrum(k,:)=A.data(l,:);
    end
end

% Remove bad rows such as focusing on an imperfection and creating 100%
% fluorescence
Spectrum(207, :) = Spectrum(206, :);
Spectrum(400, :) = Spectrum(399, :);

Spectrum(400,1)

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
    %Spectrum(i, :) -= (zeros(1, colsSpectrum) + minVal);
end

for peak_num=1:size(peaks, 1)
    for x = 0:num_steps_x-1
        for y = 1:num_steps_y
            index = (num_steps_y * x) + y;
            
            if mod(x, 2) == 0
                y_pos = num_steps_y - y + 1;
            else
                y_pos = y;
            end
            grids(x+1, y, peak_num) = sum(Spectrum(index, peaks{peak_num, 1}:peaks{peak_num, 2}));
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

for peak_num=1:size(peaks, 1)   
    subtightplot(3, 3, peak_num, [0.01, 0.03, 0.03])
    imagesc(grids(:, :, peak_num))
    pbaspect([1 num_steps_x/num_steps_y 1])
    colorbar
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

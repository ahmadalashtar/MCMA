% MCMA
function MCMAval = MCMA(Ilow, Ienh)
% Comments:
% 1- This code only works for 8-bit gray-scale images

if numel( size( Ilow ) ) == 3
    Ilow = rgb2gray( Ilow );
end
if numel( size( Ienh ) ) == 3
    Ienh = rgb2gray( Ienh );
end

[hlow, henh] = denoise(Ilow, Ienh);

PUval = getPU(Ienh);
HSDval = getHSD(hlow, henh);
DROval = getDRO( henh );
a = -.7;
b = -.3;
c = .4;

MCMAval = ((a * PUval + b * HSDval + c * DROval) + 1) / 1.4;
end
%% 



%% ---------- DRO --------------
function DROval = getDRO( henh )
idx = find(henh);
dyn = idx(end) - idx(1);
maxDyn = 256;
DROval = dyn / maxDyn;
end
%% 



%% ---------- PU -------------
function PUval = getPU(Ienh)
    Ienh = double(Ienh);

    % Calculate the size of the matrix for valid regions
    paddedIenh = padarray(Ienh, [1, 1], 'replicate');
    centerBlock = paddedIenh(2:end-1, 2:end-1);
    
    % Create a 3D block matrix for surrounding pixels
    neighbors = cell(3, 3);
    for i = 1:3
        for j = 1:3
            neighbors{i, j} = paddedIenh(i:end-3+i, j:end-3+j);
        end
    end
    
    % Stack neighbors into a single 3D matrix
    neighborStack = cat(3, neighbors{:});
    
    % Compute PU values for all pixels
    absDiff = abs(neighborStack - centerBlock);
    PU = sum(1 ./ (absDiff + 1), 3) / 9;

    % Exclude border pixels from the computation
    PU_inner = PU(2:end-1, 2:end-1);
    PUval = mean(PU_inner(:));
end

%% 




%% ---------- HSD --------------
function HSDval = getHSD(hlow, henh)

% [hlow, henh] = denoise(Ilow, Ihigh);

values_low = find(hlow);
start_low = values_low(1);
end_low = values_low(end);

values_enh = find(henh);
start_enh = values_enh(1);
end_enh = values_enh(end);

dynamic_low = hlow(start_low:end_low)';
dynamic_enh = henh(start_enh:end_enh)';

% upsample

% low 
rate=256/length(dynamic_low);
stretched_low=zeros(1,256);
for i=1:length(dynamic_low)
    destidx=ceil((i-1)*rate+1);
    stretched_low(destidx)=dynamic_low(i);
end

% interpolate
stretch_low_nz=find(stretched_low);
i_stretched_low=interp1(stretch_low_nz, stretched_low(stretch_low_nz), 1:length(stretched_low));

rate=256/length(dynamic_enh);
stretched_enh=zeros(1,256);
for i=1:length(dynamic_enh)
    destidx=ceil((i-1)*rate+1);
    stretched_enh(destidx)=dynamic_enh(i);
end

% interpolate
stretch_enh_nz=find(stretched_enh);
i_stretched_enh=interp1(stretch_enh_nz, stretched_enh(stretch_enh_nz), 1:length(stretched_enh));


nanz=find(isnan(i_stretched_low));
for i=1:length(nanz)
    i_stretched_low(nanz(i))=i_stretched_low(nanz(i)-1);
end

nanz=find(isnan(i_stretched_enh));
for i=1:length(nanz)
    i_stretched_enh(nanz(i))=i_stretched_enh(nanz(i)-1);
end

% out2 to be calculated 
% 8
stretched_low = stretched_low / sum(stretched_low);
stretched_enh = stretched_enh / sum(stretched_enh);
bin_length = 8;
diffs = zeros(bin_length,256/bin_length);
for j=1:bin_length
    for i=1:256/bin_length
        if (mod(j,bin_length)==1)
            bin_low=stretched_low((i-1)*bin_length+j:i*bin_length+j-1);
            bin_enh=stretched_enh((i-1)*bin_length+j:i*bin_length+j-1);            
        else
            if i==256/bin_length % exception
                bin_low=[stretched_low(end-bin_length+j:end), stretched_low(1:j-1)];
                bin_enh=[stretched_enh(end-bin_length+j:end), stretched_enh(1:j-1)];
            else % as usual
                bin_low=stretched_low((i-1)*bin_length+j:i*bin_length+j-1);
                bin_enh=stretched_enh((i-1)*bin_length+j:i*bin_length+j-1);
            end
        end
        diffs(j,i)=abs(sum(bin_enh)-sum(bin_low));
    end
end


HSDval = sum( sum( diffs ) ) / ( 8 * 2 ); % 8 * 2 is MAXHSD


end
%% 




%% ---------- De-noise -------------
%  --- Remove salt-pepper noinse----
function [denoisedHLow , denoisedHEnh] = denoise( Ilow , Ienh )

hlow = imhist(uint8(Ilow));
henh = imhist(uint8(Ienh));

% remove salt & pepper - threshold = 0.001
thr = 0.001;
pix_thr = round(numel(Ilow) * thr);
% remove salt 
idx = 1;
for i = 2 : 256
    if sum(henh(end - i + 1 : end)) < pix_thr
        idx = i;
    else
        break;
    end
end
if idx > 1
    henh(end - idx + 1 : end) = 0;
end

idx = 1;
for i = 2 : 256
    if sum(hlow(end - i + 1 : end)) < pix_thr
        idx = i;
    else
        break;
    end
end
if idx > 1
    hlow(end - idx + 1 : end) = 0;
end

% remove pepper
idx = 1;
for i = 2 : 256
    if sum(henh(1 : i)) < pix_thr
        idx = i;
    else
        break
    end
end
if idx > 1
    henh(1 : idx) = 0;
end

idx = 1;
for i = 2 : 256
    if sum(hlow(1 : i)) < pix_thr
        idx = i;
    else
        break;
    end
end
if idx > 1
    hlow(1 : idx) = 0;
end

denoisedHLow = hlow;
denoisedHEnh = henh;

end
%%----------------------------------
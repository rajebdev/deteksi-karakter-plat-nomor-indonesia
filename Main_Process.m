close all;
clear;
clc;

gambar = imread('plat2.jpg');

ukuran = size(gambar);

%================ PREPROCESSING ============

%Resize Gambar
imgSmall = imresize(gambar, 0.6);

%convert to grayscale
imgGray = rgb2gray(imgSmall);

%Deteksi Tepi dengan menggunakan sobel
imgEdge = edge(imgGray, 'sobel');

%Menggandeng Pixel yang kurang dari 10
s10 = strel('disk', 10);
imgCloseEdge = imclose(imgEdge, s10);

% Menampilkan Gambar Biar Good
figure;
subplot(2,1,1);
imshow(imgSmall); title ('Gambar Asli');
subplot(2,1,2);
imshow(imgGray); title('Gambar Abu Abu');

% Menghilangkan Noise (10 pixel)
imgClearNoise = bwareaopen(imgCloseEdge, 10);

%Menggandeng pixel yang letaknya kurang dari 2 pixel
s2 = strel('disk', 2);
imgClosePixel = imclose(imgClearNoise, s2);

%Menghilangkan Lubang area yang tertutup
imgClearHole = imfill(imgClosePixel, 'holes');

%Mencari Garis Pinggir
[B, L] = bwboundaries(imgClearHole, 'noholes');
maxArea = 0;

stats = regionprops(L, 'Area', 'Centroid', 'BoundingBox', 'Orientation');

threshold = 0.90;

%Looping Boundaries
for k=1:length(B)
    
    %Mencari B (X,Y) label pada index K
    boundary = B{k};

    %Menghitung Keliling Object
    delta_sq = diff(boundary).^2;
    perimeter = sum(sqrt(sum(delta_sq, 2)));
    
    %Mencari Luas Area Object
    area = stats(k).Area;
    
    %Menghitung metric seberapa dekat dengan bentuk lingkaran
    metric = 4*pi*area/perimeter^2;
    
    %Menampilkan Hasil Perhitungan
    areatulis  = fprintf('Area %2.2f\n', area);
    if area > maxArea
        maxArea = area;
    end
    
    %Menandai Centroid
    if metric > threshold
        centroid = stats(k).Centroid;
    end
end

fprintf('Max %2.2f\n', maxArea);

%Mencari Luas Area Terluas
for k = 1: length(B)
    area = stats(k).Area;
    if maxArea == area
        cropResult = stats(k).BoundingBox;
    end
end

%Croping Gambar
imgOri = imcrop(imgSmall, cropResult);
imgGray = imcrop(imgGray, cropResult);
imgEdge = imcrop(imgEdge, cropResult);
imgClear = imcrop(imgClearHole, cropResult);


imgClear2 = imresize(imgOri, [600 1000], 'bilinear');
imgGray = imresize(imgGray, [600 1000], 'bilinear');
%imgThumb = imcrop(imgGray,[250 100 200 180]);
imgThumb = imbinarize(imgGray, 0.65);

figure;
imshow(imgThumb), title('Hasil Croping');

imgThumb = reshape(imgThumb, 1, []);
imgThumb = mean(imgThumb);

%Mendetekti apabila gambar terlalu kecil
if imgThumb >= 0.65
    imgGray = imcomplement(imgGray);
    figure;
    imshow(imgGray);
end


%Mencari Warna plat nmor
imgAdaptive = adaptthresh(imgGray, 25*0.01);
imgAdaptive = imbinarize(imgGray, imgAdaptive);
imgWarnaBalik = imcomplement(imgAdaptive);
imgBw = bwareaopen(imgWarnaBalik, 2000);

figure, imshow(imgBw), title('Gambar Warna Balik');

%Mencari garis pinggir%
[BB,LL] = bwboundaries(imgBw,'holes');
maxArea=0

stats = regionprops(LL,'Area','Centroid','BoundingBox','Orientation');

% loop boundaries
for kk = 1:length(BB)
    % Mencari koordinat (X,Y)label 'k'
    boundary = BB{kk};
    
    % object's perimeter
    delta_sq = diff(boundary).^2;
    perimeter = sum(sqrt(sum(delta_sq,2)));
    
    % Mencari luas area
    area = stats(kk).Area;
    % Menghitung roundness metric
    metric = 4*pi*area/perimeter^2;
    
    % Menampilkan hasil
    % metric_string = sprintf('%2.2f',metric);
    fprintf('%2.2f',area);
    
    if area > maxArea
        maxArea=area;
    end
    % Menandai
    if metric > threshold
        centroid = stats(kk).Centroid;
    end
end


%Membuat Matrik Baru Hasil Normalisasi
Matrik = [];
n=1;
figure
for kk = 1:length(BB)
    area = stats(kk).Area;
    if area >= 1000 && area <= 23000
        cropResult2 = stats(kk).BoundingBox;
        imgResult = imcrop(imgBw, cropResult2);
        [tinggi, panjang] = size(imgResult);
        fprintf('Tinggi = %2.2f\nPanjang = %2.2f\n', tinggi, panjang);
        rasio = tinggi/panjang;
        fprintf('Rasio %2.2f\n', rasio);
        if tinggi > 130 && rasio >= 0.85 && rasio <= 6
            subplot(1,9,n);imshow(imgResult);title(area);
            Matrik = [Matrik tinggi];
            n=n+1;
        end
    end
end
medianMatrik = median(Matrik);
stdMatrik = std(Matrik);

%Mencari Perbandingan Tinggi dan Panjang
figure
n=1;
for kk =1:length(BB)
    area = stats(kk).Area;
    if area >= 1000 && area <= 23000
        cropResult3 = stats(kk).BoundingBox;
        imgResult = imcrop(imgBw, cropResult3);
        [tinggi, panjang] = size(imgResult);
        if tinggi >= 130
            rasio = tinggi/panjang;
            if rasio >= 0.5 && rasio <= 10
                if tinggi >= medianMatrik - (stdMatrik+10)
                    if tinggi <= medianMatrik + (stdMatrik+10)
                        imgResult = imclose(imgResult, strel('diamond', 1));
                        imgResult = imresize(imgResult, [20 10], 'bilinear');
                        imgResultRe = im2bw(imgResult, 0.4);
                        subplot(1, 9, n); imshow(imgResultRe), title(area);
                        n=n+1;
                    end
                end
            end
        end
    end
end

                        
                        










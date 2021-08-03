clc; clear; close all;

tic

%% Inputs

rng('default') %for reproducability



%% Background

MaskImage = double(imread('Mask.tif'))/(2^8-1);
x1=find(sum(MaskImage,2),1);
x2=find(sum(MaskImage,2),1,'last');
y1=find(sum(MaskImage,1),1);
y2=find(sum(MaskImage,1),1,'last');





%% The bulk

for N = 5:5:30  %pixels/pixel
    BackgroundImage = double(imread(strcat('MAX_C3-segmentation_002_C.tif')))/(2^16-1).*MaskImage;
    BackgroundImage = BackgroundImage(x1:x2,y1:y2,:);
    BackgroundImage = scaleto(BackgroundImage,N);
    for DontIncludeCell = [true,false]
        Folder = 'gene_nucleus/';
        File = dir(Folder);
        numGenePics = length(File)-3;
        Gene = {};
        for ge = 1:numGenePics
            file = strcat(Folder,File(ge+3).name);

            I = double(imread(file))/(2^16-1).*MaskImage;
            I = I(x1:x2,y1:y2,:);
        %     figure; imshow(I)
            I = scaleto(I,N);
            Gene{ge} = I;
        end


        if ~DontIncludeCell
            numGenePics = numGenePics + 1;
            Gene{numGenePics} = BackgroundImage;
        end



        genepics = zeros(numel(Gene{1}),numGenePics);
        for ge = 1:numGenePics
            genepics(:,ge) = Gene{ge}(:);
        end

        badrows = all(genepics == 0,2);
        goodrows = ~all(genepics == 0,2);
        genepics(badrows,:) = [];

        opts = statset('MaxIter',10^3);
        [Y, loss] = tsne(genepics,'Distance','correlation','NumDimensions',2,'LearnRate',5000);%,'Options',opts,'NumPrint',1,'Verbose',1);
        Z = linkage(Y,'ward');
        for colorsinpic = [3:10]
            c = cluster(Z,'Maxclust',colorsinpic);
            % scatter(Y(:,1),Y(:,2),10,c)
            cnew = zeros(numel(Gene{1}),1);
            cnew(goodrows) = c;
            newim = reshape(cnew,size(Gene{1}));

            I = label2rgb(newim, 'hsv','k','shuffle');

            % figure; imshow(I)

            NamePart = 'Without';
            TSNEfile = strcat('TSNEpic',NamePart(1:4 + DontIncludeCell*3),'Cell',num2str(colorsinpic),'Pic',num2str(N),'.jpg');
            imwrite(I,TSNEfile)
            
            
            
            
            CountTheBrightness = zeros(colorsinpic,size(genepics,2));
            for ii = 1:colorsinpic
                CountTheBrightness(ii,:) = mean(genepics.*(cnew(goodrows) == ii),1);
            end
            ClusterNumberHere = [1:colorsinpic].';
            GeneNumberHere = [1:size(genepics,2)].';
            
            
            GENES = {};
            for ii = 1:length(GeneNumberHere)
                if ii <= 7
                    nameofgenehere = File(ii+3).name;
                    throwaway = find(nameofgenehere == '_');
                    throwaway = throwaway(end)+1;
                    aaaaab = find(nameofgenehere == '.')-1;
                    nameofgenehere = nameofgenehere(throwaway:aaaaab);
                    GENES{ii} = nameofgenehere;
                else
                    GENES{ii} = 'Cell';
                end
            end
            h = heatmap(GENES,ClusterNumberHere,CountTheBrightness);
            h.Title = 'Brightness Of Genes By Color';
            h.XLabel = 'Cluster';
            HEATMAPfile = strcat('TSNEclusterHeatmap',NamePart(1:4 + DontIncludeCell*3),'Cell',num2str(colorsinpic),'Pic',num2str(N),'.jpg');
            saveas(gcf,HEATMAPfile);

        end
        
        
        
        
    end
end




%% Count Dots

Folder = 'gene_nucleus/';
File = dir(Folder);
numGenePics = length(File)-3;
Gene = {}; geneName = {};
for ge = 1:numGenePics
    
    nameofgenehere = File(ge+3).name;
    throwaway = find(nameofgenehere == '_');
    throwaway = throwaway(end)+1;
    aaaaab = find(nameofgenehere == '.')-1;
    geneName{ge} = nameofgenehere(throwaway:aaaaab);
    
    file = strcat(Folder,File(ge+3).name);

    I = double(imread(file))/(2^16-1).*MaskImage;
    I = I(x1:x2,y1:y2,:);
    B = imgaussfilt(I,1/2) - imgaussfilt(I,2);
    B = B.*(B>0);
    J = (B>multithresh(B));
    
    
    S = sum(bwmorph(J,'shrink',inf),'all');
    Gene{ge} = S;
end




%% We're done!!!
disp('Done!!!')

RuntimeSeconds = toc
RuntimeMinutes = toc/60
RuntimeHours = toc/60^2




%% Functions used!!!
function [I] = scaleto(I,N)
    [rows, columns, ~] = size(I);
    numOutputRows = floor(rows/N);
    numOutputColumns = floor(columns/N);
    I = imresize(I(1:N*floor(end/N),1:N*floor(end/N),:), [numOutputRows, numOutputColumns]);

end






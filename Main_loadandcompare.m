% This file directly compares results from mol
clear all, close all

%% load data
rootpath = 'C:\Users\LinR\OneDrive - University of Twente\work\project\Spectral analysis of HRV, PTT and BP\Data\AST\Data';
fileInfo = {};
nn=0; count=0;
dirInfo = dir(rootpath);
young_idx = [60 61 63 64];
old_idx = [36 38 39 40 42 45 49];
for i = 1:length(dirInfo)
    if ~dirInfo(i).isdir || startsWith(dirInfo(i).name, '.')
        continue;
    end

    firstLevelFolder = dirInfo(i).name;
    firstLevelPath = fullfile(rootpath, firstLevelFolder);

    subDirInfo = dir(firstLevelPath);

    for j = 1:length(subDirInfo)
        if ~subDirInfo(j).isdir || startsWith(subDirInfo(j).name, '.')
            continue;
        end

        secondLevelPath = fullfile(firstLevelPath, subDirInfo(j).name);

        thirdDirInfo = dir(secondLevelPath);

        for k = 1:length(thirdDirInfo)
            if ~thirdDirInfo(k).isdir || startsWith(thirdDirInfo(k).name, '.')
                continue;
            end

            thirdLevelFolder = thirdDirInfo(k).name;
            thirdLevelPath = fullfile(secondLevelPath, thirdLevelFolder);

            derivedFileName = sprintf('%s_%s_preprocessed.mat', firstLevelFolder, thirdLevelFolder);
            derivedFilePath = fullfile(thirdLevelPath, derivedFileName);

            annotationFileName = sprintf('%s_annotations.mat', firstLevelFolder);
            annotationFilePath = fullfile(thirdLevelPath, annotationFileName);
            
            load([thirdLevelPath '\' annotationFileName])
            
            if isfield(annotation,"HGS")  %|| isfield(annotation,"S_to_stand") 
                % disp(derivedFileName)
                % nn=nn+1;
            load([thirdLevelPath '\'  derivedFileName])
            [output,~,count_]=processing_posturalchange(data,annotation,dirInfo(i).name);
            end
            % count  = count_ + count;
            % if isfield(output,'HR')
            %     HR_feature.SDNN = [HR_feature.SDNN output.HR.SDNN]; HR_feature.SDSD = [HR_feature.SDSD output.HR.SDSD]; 
            %     HR_feature.pNN50 = [HR_feature.pNN50 output.HR.pNN50]; HR_feature.RMSSD = [HR_feature.RMSSD output.HR.RMSSD];             
            %     HR_feature.ApEn = [HR_feature.ApEn output.HR.ApEn]; HR_feature.triangular_val = [HR_feature.triangular_val output.HR.triangular_val];
            %     HR_feature.SD1 = [HR_feature.SD1 output.HR.SD1]; HR_feature.SD2 = [HR_feature.SD2 output.HR.SD2];
            %     HR_feature.pLF = [HR_feature.pLF output.HR.pLF]; HR_feature.triangular_val = [HR_feature.pHF output.HR.pHF];
            %     HR_feature.LFHFratio = [HR_feature.LFHFratio output.HR.LFHFratio_HR]; HR_feature.SDratio = [HR_feature.SDratio output.HR.SDratio];
            % 
            %     PR_feature.SDNN = [PR_feature.SDNN output.PR.SDNN]; PR_feature.SDSD = [PR_feature.SDSD output.PR.SDSD]; 
            %     PR_feature.pNN50 = [PR_feature.pNN50 output.PR.pNN50]; PR_feature.RMSSD = [PR_feature.RMSSD output.PR.RMSSD];             
            %     PR_feature.ApEn = [PR_feature.ApEn output.PR.ApEn]; PR_feature.triangular_val = [PR_feature.triangular_val output.PR.triangular_val];
            %     PR_feature.pLF = [PR_feature.pLF output.PR.pLF]; PR_feature.pHF = [PR_feature.pHF output.PR.pHF];
            %     PR_feature.SD1 = [PR_feature.SD1 output.PR.SD1]; PR_feature.SD2 = [PR_feature.SD2 output.PR.SD2];
            %     PR_feature.LFHFratio_PR = [PR_feature.LFHFratio output.PR.LFHFratio_PR]; PR_feature.SDratio = [PR_feature.SDratio output.PR.SDratio];
            % end
            % save('Results.mat',"PR_feature","HR_feature")
            % else, continue; end
        end
    end
end
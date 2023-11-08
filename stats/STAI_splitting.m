%%%
%   NeuroTin / 20230116 / STAI_splitting.m
%
%   Split the STAI questionnaire answers (from Evamed
%   — the online reporting platform) into STAI Y-1 and
%   STAI Y-2 data for all participants (CBT, fMRI NF,
%   and EEG NF arms).
%%%
%
%   Nicolas Gninenko / nicolas.gninenko@gmail.com
%
%%%%


% Load data
%evamedSTAI = readcell('Evamed_wyss-neurotin__2023-01-14.xlsx');
evamedSTAI = readcell('Evamed_wyss-neurotin__2023-01-14_Evamed_int.xlsx');

CBT_patients = load('CBT_fMRI_patients_cells.mat','CBT_patients');
CBT_patients = CBT_patients.CBT_patients;
fMRI_patients = load('CBT_fMRI_patients_cells.mat','fMRI_patients');
fMRI_patients = fMRI_patients.fMRI_patients;


% Fill in the data for STAI separately
STAI_Y1_CBT = zeros(length(CBT_patients),3); % column 1: baseline; 2: post early (+1 mo.); 3: post late (+6 mo.)
STAI_Y2_CBT = zeros(length(CBT_patients),3);
STAI_Y1_fMRI = zeros(length(fMRI_patients),3);
STAI_Y2_fMRI = zeros(length(fMRI_patients),3);

cc = 1;
for p = fMRI_patients'
    tmp_line = find(cell2mat(cellfun(...
        @isequal,evamedSTAI(:,4),repmat(p,size(evamedSTAI,1),1),'UniformOutput',false)),Inf);
    if length(tmp_line)~=3 && ~strcmp(p{1},'025GNV'), error('Check manually.'); end % handle special case
    for col = 14:33
        if ~isnumeric(evamedSTAI{tmp_line(1),col})
            STAI_Y1_fMRI(cc,1) = NaN; break;
        else
            STAI_Y1_fMRI(cc,1) = STAI_Y1_fMRI(cc,1) + evamedSTAI{tmp_line(1),col};
        end
    end
    for col = 34:53
        if ~isnumeric(evamedSTAI{tmp_line(1),col})
            STAI_Y2_fMRI(cc,1) = NaN; break;
        else
            STAI_Y2_fMRI(cc,1) = STAI_Y2_fMRI(cc,1) + evamedSTAI{tmp_line(1),col};
        end
    end
    for col = 60:79
        if ~isnumeric(evamedSTAI{tmp_line(2),col})
            STAI_Y1_fMRI(cc,2) = NaN; %break;
        else
            STAI_Y1_fMRI(cc,2) = STAI_Y1_fMRI(cc,2) + evamedSTAI{tmp_line(2),col};
        end
        if length(tmp_line)==3
            if ~isnumeric(evamedSTAI{tmp_line(3),col})
                STAI_Y1_fMRI(cc,3) = NaN; %break;
            else
                STAI_Y1_fMRI(cc,3) = STAI_Y1_fMRI(cc,3) + evamedSTAI{tmp_line(3),col};
            end
        else
            STAI_Y1_fMRI(cc,3) = NaN; %break;
        end
    end
    for col = 80:99
        if ~isnumeric(evamedSTAI{tmp_line(2),col})
            STAI_Y2_fMRI(cc,2) = NaN; %break;
        else
            STAI_Y2_fMRI(cc,2) = STAI_Y2_fMRI(cc,2) + evamedSTAI{tmp_line(2),col};
        end
        if length(tmp_line)==3
            if ~isnumeric(evamedSTAI{tmp_line(3),col})
                STAI_Y2_fMRI(cc,3) = NaN; %break;
            else
                STAI_Y2_fMRI(cc,3) = STAI_Y2_fMRI(cc,3) + evamedSTAI{tmp_line(3),col};
            end
        else
            STAI_Y2_fMRI(cc,3) = NaN; %break;
        end
    end
    cc = cc+1;
end

cc = 1;
for p = CBT_patients'
    tmp_line = find(cell2mat(cellfun(...
        @isequal,evamedSTAI(:,4),repmat(p,size(evamedSTAI,1),1),'UniformOutput',false)),Inf);
    if length(tmp_line)~=3 && ~strcmp(p{1},'064GNV'), error('Check manually.'); end % handle special case
    for col = 14:33
        if ~isnumeric(evamedSTAI{tmp_line(1),col})
            STAI_Y1_CBT(cc,1) = NaN; break;
        else
            STAI_Y1_CBT(cc,1) = STAI_Y1_CBT(cc,1) + evamedSTAI{tmp_line(1),col};
        end
    end
    for col = 34:53
        if ~isnumeric(evamedSTAI{tmp_line(1),col})
            STAI_Y2_CBT(cc,1) = NaN; break;
        else
            STAI_Y2_CBT(cc,1) = STAI_Y2_CBT(cc,1) + evamedSTAI{tmp_line(1),col};
        end
    end
    for col = 60:79
        if ~isnumeric(evamedSTAI{tmp_line(2),col})
            STAI_Y1_CBT(cc,2) = NaN; %break;
        else
            STAI_Y1_CBT(cc,2) = STAI_Y1_CBT(cc,2) + evamedSTAI{tmp_line(2),col};
        end
        if length(tmp_line)==3
            if ~isnumeric(evamedSTAI{tmp_line(3),col})
                STAI_Y1_CBT(cc,3) = NaN; %break;
            else
                STAI_Y1_CBT(cc,3) = STAI_Y1_CBT(cc,3) + evamedSTAI{tmp_line(3),col};
            end
        else
            STAI_Y1_CBT(cc,3) = NaN; %break;
        end
    end
    for col = 80:99
        if ~isnumeric(evamedSTAI{tmp_line(2),col})
            STAI_Y2_CBT(cc,2) = NaN; %break;
        else
            STAI_Y2_CBT(cc,2) = STAI_Y2_CBT(cc,2) + evamedSTAI{tmp_line(2),col};
        end
        if length(tmp_line)==3
            if ~isnumeric(evamedSTAI{tmp_line(3),col})
                STAI_Y2_CBT(cc,3) = NaN; %break;
            else
                STAI_Y2_CBT(cc,3) = STAI_Y2_CBT(cc,3) + evamedSTAI{tmp_line(3),col};
            end
        else
            STAI_Y2_CBT(cc,3) = NaN; %break;
        end
    end
    cc = cc+1;
end


% Some statistics (temp)
% [p_,~] = vartestn([[STAI_Y1_fMRI;[]] STAI_Y1_CBT],...
%     'TestType','LeveneQuadratic'); % variance test across subgroups for Y1
% fprintf(['\nSTAI-Y1 fMRI vs CBT (PA, EA, LA): ' num2str(p_) '.']);
% [p_,~] = vartestn([[STAI_Y2_fMRI;[]] STAI_Y2_CBT],...
%     'TestType','LeveneQuadratic'); % variance test across subgroups for Y2
% fprintf(['\nSTAI-Y2 fMRI vs CBT (PA, EA, LA): ' num2str(p_) '.\n']);
%close all;




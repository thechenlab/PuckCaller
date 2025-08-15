
%%%%This pipeline function assumes that we are using the large image
%%%%feature and the XY feature in Nikon images, and that files ending in
%%%%xy1 are for puck 1. We do not do stitching.

%This function is for beads with 7 J bases in two groups.

%%%%SETUP:

%0) Make sure ImageJ is closed.

%1) When you are done with the code, you should make the raw data folder in
%Pucks and the InputFolder in find_roi online only via smart sync to save
%hard drive space. You should also delete the pucktmp directory.

%2) A number of paths are hardcoded in the code currently, e.g.
%"C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers\vlfeat-0.9.20\toolbox\"
%in find_roi_stack_fun

%3) Change ImageSize to reflect something slightly smaller than the size of
%the final stitched images.

%4) If you have a digital gene expression matrix, you need to take the entire CSV and put it in the Illumina folder (C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\Puck_180106_3) and call it DGE.csv. I then copy out the top row with notepad++, delete the first entry (which is there as a placeholder for the row labels) and put it into a different file called IlluminaBarcodes.txt

%5) Each ligation sequence and primer set will require a bunch of modifications
%to MapLocationsFun. To accommodate for this, I make a different version of MapLocationsFun for
%each ligation sequence. You have to change which version is called below.

%6) NOTE: SlideseqSave.ijm is stored in C:\Fiji.app\macros, and the text
%is:
%params=getArgument()
%print("Opening file:")
%print(params)
%open(params);
%saveAs("Tiff","C:\\PuckOutputTmp.tif")
%print("Done.")
%exit();


%For future releases, consider setting the priority of the parpool
%processes to High:

%matlabpool open
%cmd_str = 'wmic process where name="MATLAB.exe" CALL setpriority 64';
%[~,~] = system(cmd_str);


%% Initialize
%Set up the paths to all of the various folders with functions in them that
%we will call:
clear all
close all

PythonPath='C:\Users\sgr\AppData\Local\Programs\Python\Python37\python.exe';
BeadseqCodePath='\\helium\broad_thechenlab\Gio\BeadSeq Code';
PipelineFunctionPath='\\helium\broad_thechenlab\Gio\PipelineFunctions';
addpath(BeadseqCodePath,[BeadseqCodePath,'\find_roi'],PipelineFunctionPath);
addpath([BeadseqCodePath,'\find_roi\helpers']);
%We assume that the nd2s have been exported to tiffs in the format:
%DescriptiveNametLYX, where DescriptiveName refers to one run of the microscope,  Y is a letter and X is a number, and L is the
%name of the ligation within that DescriptiveName file series.
%We convert the files to a final output format:
%Puck85 Ligation X Position AB
%BeadType="180402"; for 14bp barcodes from 180402 beads
BeadType="180728";

RenameFiles=1;

RunSignificanceAnalysis=1;  %not important
NumClusters=[11,17,14,14]; %Cerebellum is 1, hippocampus is 2, frontalcortex is 3, posteriorcortex is 4 % not important 

%Which pucks do we try to run analogizer on? 0 if don't run, otherwise, the
%number here refers to the element of DropseqDGEPaths and
%DropseqClusterPaths to use as a reference.
RunAnalogizer=[0,2,2,0,2,0,2,2,2,2,2,2,2,0,0,0,1,1,1,0,2,2,2,2]; %not i,portant
AnalogizerBeadCutoff=5; %Note that this threshold is only applied for the variable genes used in the NMFreg, so it can be low
AnalogizerType="NMFReg";
%Cerebellum is 1. Hippocampus is 2.
DropseqDGEPaths={'\\helium\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT','\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Hippocampus'};
DropseqClusterPaths={'\\helium\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT\assign\F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS','\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Hippocampus\assign\F_GRCm38.81.P60Hippocampus.cluster.assign.RDS'};
DropseqMeanAndVariancePath=fullfile(BeadseqCodePath,'DGEMeansAndVariances');

NumPar=20; %number of threads

CropImage=1; %flag that asks you to crop output. THis is not used 
CropSuffix='_Cropped';

EnforceBaseBalance=1; 
BaseBalanceTolerance=0.05;

SaveData=1;
IlluminaReadThreshold=10;
MaximumBarcodesToAnalyze=160000;
NumLigations=20; %We assume that the missing ligation is Truseq-4 if you are doing 13.
NumBases=14; %This is the number of bases sequenced
%BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14];
%BarcodeSequence=[1,2,3,4,0,5,0,6,7,8,9,10,0,11,0,12,0,13]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
tmpfolder='\\helium\broad_thechenlab\Gio\pucktmp\';

ImageSize=[8400, 8400]; %The image registration will output images that are not all the exact same size, because of stitching. So find_roi_stack_fun crops the images a bit. We can choose this value to be something like 0.95 * the size of the images. So e.g. for 3x3 it should be 0.95*2048*(2*(2/3)+1) = 4500. For 7x7 it is 0.95*2048*(1+6*2/3)=9728. I used 10400 previously though.
XCorrBounds=[4000,4250,4000,4250]; %This is the ROI used in channel registration
RegisterColorChannels=1;
BeadZeroThreshold=1;
PixelCutoffRegistration=400;
PixelCutoffBasecalling=300;
DropBases=0;
BeadSizeCutoff=7; %changed from 30 to 8 for downsample

%The illumina barcodes are assumed to be 13 bases long, consisting of the
%first 6 J bases and then the last 7 J bases. If the barcodes you are using
%are different, you have to change it in the bs2cs calls.

%And ligation sequnece:
%Tru_L1, Tru_L2, Tru-1_L1, Tru-1_L2, Tru-2_L1, Tru-2_L2, Tru-3_L1, Tru-3_L2, Tru-4_L1, Tru-4_L2
%The numbers in these variables are the ***second base being interrogated***
%Equivalently, it is the number associated with that ligation in figure 4a of the Fisseq nature protocols paper
PrimerNLigationSequence = [2, 7, 1, 6, 5, 4, 3]; %good for 14 ligations
%PrimerNLigationSequence = [2, 7, 1, 6, 5, 4]; %good for 13 ligations
%UP_L1, UP_L2, UP-1_L1, UP-1_L2, UP-2_L1, UP-2_L2, UP-3_L1, UP-3_L2, UP-4_L1, UP-4_L2
PrimerUPLigationSequence=[2, 7, 1, 6, 5, 4,3];

%Note that bs2cs will output the colors corresponding to each ligation in this ligation sequence, in the order specified here,
%So if you were only to do 6 ligations off primer N instead of 7, you would only need to
%remove the final element in the PrimerNLigationSequence

InverseLigationSequence=[3,1,7,6,5,4,2,10,8,14,13,12,11,9]; %Good for both 13 and 14 ligations.
%Before exporting, "N"s are added to the end of the ligation
%sequence, and must then be inserted into the correct location to stand in
%for the missing ligations. This code puts the N at position 7
%WhichLigationsAreMissing=[1,2,3,4,5,6,14,7,8,9,10,11,12,13];
WhichLigationsAreMissing=[1,2,3,4,5,6,7,8,9,10,11,12,13,14];

%Mask
Mask = logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
%ogical mask for which ligation bases are bad and should not be used during
%mapping

%We need to determine if the device running the pipeline must be linked to
%the google doc for updating basecalling metrics into the puck log. Note
%that this must be set up for each unique device that runs the pipeline.
%The user will be prompted to manually enter a given security code into
%google.com/devices. After this is done once for a given device, subsequent
%pipeline runs can update automatically without any user action. Of note,
%privacy setting on puck log document must be set to "anyone with link can
%edit".

% reply = input('Is this the first time you will be running the Slide-seq pipeline on this device since 3/25/19 (reply: y/n)?','s');
% if strcmp(reply,'y')
%   %RunOnce(clientID,client_secret)
%   RunOnce('248080687373-vc0mevs6chk2rgpm44uva797tki6upf4.apps.googleusercontent.com', '6uYApyl6mlkwIk3cNH_qUVl')
% end


%% Load parameters from manifest file
display('Load parameters from manifest file');

% Input/select manifest file such as 'D:\Jilong\out\manifest.v2.txt'
% prompt = 'Please input the full path of your manifest file:\n';
% manifest = input(prompt);
[file,path] = uigetfile('*.txt', 'Select manifest file');
manifest = [path,file];

% Check if manifest file exists
if ~exist(manifest,'file')
    error('Manifest file not found');
end

% Load manifest content into a variable set
fid = fopen(manifest, 'rt');
fcon = textscan(fid,  '%s%s', 'Delimiter', '=');
fclose(fid);
indat = horzcat(fcon{:});
[m,n] = size(indat);
params = [];
for i=1:m
   params.(indat{i,1}) = indat{i,2}; 
end

% Retrieve variable values based on variable names
Nd2Folder = params.('Nd2Folder');
FolderWithProcessedTiffs = params.('FolderWithProcessedTiffs');
OutputFolderRoot = params.('OutputFolderRoot');
IndexFiles = params.('IndexFiles');
PuckName = params.('PuckName'); % Puck_190327
PucksToAnalyze = params.('PucksToAnalyze');
LigationToIndexFileMapping = params.('LigationToIndexFileMapping');
tnumMapping = params.('TnumMapping');
DeleteIntermediateFiles = 'false';
if isfield(params,'DeleteIntermediateFiles')
    DeleteIntermediateFiles = params.('DeleteIntermediateFiles');
end

PuckImageSubstraction = 'True';
if isfield(params,'PuckImageSubstraction')
    PuckImageSubstraction = params.('PuckImageSubstraction');
end

Monobase = 0;

if isfield(params,'Monobase')
    Monobase = textscan(params.('Monobase'),'%f','Delimiter',',');
    Monobase = transpose(horzcat(Monobase{:}));
    Monobase = logical(Monobase);
end

if Monobase==1
    BaseBalanceTolerance=0.01;
end

% Convert string to cell array and vector
IndexFiles = textscan(IndexFiles,'%s','Delimiter',',');
IndexFiles = horzcat(IndexFiles{:});
LigationToIndexFileMapping = textscan(LigationToIndexFileMapping,'%f','Delimiter',',');
LigationToIndexFileMapping = horzcat(LigationToIndexFileMapping{:});
tnumMapping = textscan(tnumMapping,'%f','Delimiter',',');
tnumMapping = horzcat(tnumMapping{:});
PucksToAnalyze = textscan(PucksToAnalyze,'%f','Delimiter',',');
PucksToAnalyze = horzcat(PucksToAnalyze{:});

BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14];
if isfield(params,'MissingBarcodeSequence')
    MissingBarcodeSequence = params.('MissingBarcodeSequence');
    MissingBarcodeSequence = textscan(MissingBarcodeSequence,'%f','Delimiter',',');
    MissingBarcodeSequence = horzcat(MissingBarcodeSequence{:});
    BarcodeSequence=[];
    j = 1;
    for i=1:20
        if ismember(i, MissingBarcodeSequence)
            BarcodeSequence=[BarcodeSequence;0];
        else
            BarcodeSequence=[BarcodeSequence;j];
            j = j + 1;
        end
    end
end

% Create TnumMapping from LigationToIndexFileMapping
% LigationToIndexFileMapping=[1,2,3,3,3,3,3,3,4,5,6,6,6,6,6,6,6,6,6,6];
% tnumMapping=[1,1,1,2,3,4,5,6,1,1,1,2,3,4,5,6,7,8,9,10];
%m=size(LigationToIndexFileMapping);
%tnumMapping=LigationToIndexFileMapping;
%k=1;
%for i=2:m
%    if LigationToIndexFileMapping(i)~=LigationToIndexFileMapping(i-1)
%        k=0;
%    end
%    k=k+1;
%    tnumMapping(i)=k;
%end

% Create PuckNames from PuckName and PuckSToAnalyze
% Note that the order in PuckNames should match the order in the .nd2 file.
[m,n]=size(PucksToAnalyze);
PuckNames=string(m);
for i=1:m
    PuckNames(i)=[PuckName,'_',pad(num2str(PucksToAnalyze(i)),2,'left','0')];
end


%% Create folders
OutputFolders={};
for puck=1:length(PuckNames)
    ProcessedImageFolders{puck}=[FolderWithProcessedTiffs,PuckNames{puck},'\'];
    mkdir([FolderWithProcessedTiffs,PuckNames{puck}]);
    OutputFolders{puck}=[OutputFolderRoot,PuckNames{puck},'\'];
    mkdir(OutputFolders{puck});
    
    % Copy manifest file to output folder
    copyfile(manifest,[OutputFolders{puck},'\']);
end


%% Convert .nd2 to .tif
display('Run bfconvert to convert .nd2 to .tif');

for pucknum=1:length(PuckNames)
    puck=PucksToAnalyze(pucknum);
    for ligation=1:NumLigations
        if BarcodeSequence(ligation)==0
            continue;
        end
        
        tnum=tnumMapping(ligation);
        filename=[Nd2Folder,IndexFiles{LigationToIndexFileMapping(ligation)},'.nd2'];
        if ~exist(filename,'file')
            display(['file',32,filename,32,'not found']);
            continue;
        end
        
        % Run showinf to check if puck and tnum are valid
        r = randi([10000000 99999999],1);
        outputfilename=[tmpfolder,'showinf_',num2str(r),'.txt'];
        commandfile=fopen('C:\showinfCommand.cmd','w');
        fwrite(commandfile,strcat('C:\bftools\showinf',32,'"',filename,'"',32,'>',32,outputfilename));
        fclose(commandfile);
        !C:/showinfCommand
        text = fileread(outputfilename);
        text = regexp(text, '\n', 'split');
        IdxSeriesCount = find(contains(text,'Series count ='));
        IdxSizeT = find(contains(text,'SizeT ='));
        IdxSizeT = IdxSizeT(1);
        SeriesCount = split(text(IdxSeriesCount),'=');
        SizeT= split(text(IdxSizeT),'=');
        SeriesCount = strtrim(SeriesCount(2));
        SizeT= strtrim(SizeT(2));
        SeriesCount = str2num(SeriesCount{1});
        SizeT = str2num(SizeT{1});
        
        if puck > SeriesCount
            display('puck is invalid');
            continue;
        end
        if tnum > SizeT
            display('tnum is invalid');
            continue;
        end
        delete(outputfilename);
        
        % Convert .nd2 to .tiff
        outputfilename=[ProcessedImageFolders{pucknum},PuckNames{pucknum},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.tif'];
        commandfile=fopen('C:\bfconvertCommand.cmd','w');
        % Java heap size issue without using '-tilex 512 -tiley 512'
        fwrite(commandfile,strcat('C:\bftools\bfconvert -tilex 512 -tiley 512 -series',32,num2str(puck-1),' -timepoint',32,num2str(tnum-1),' "',replace(filename,'\','\\'),'" "',replace(outputfilename,'\','\\'),'"'));
        fclose(commandfile);
        !C:/bfconvertCommand

        % Convert .tiff to .png
        %outputfilename2=[ProcessedImageFolders{pucknum},PuckNames{pucknum},'_Ligation_',pad(num2str(ligation),2,'left','0'),'_Stitched.png'];
        %tiffile = imread(outputfilename, 2);
        %imwrite(tiffile, outputfilename2);
        %delete(outputfilename);
    end
end


%% Registration

display('Image Registration')

for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
    display(['Beginning registration on puck number ',num2str(puck)])
    %We now send the command to do the registration with find_roi
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    Suffix='_Stitched';
    find_roi_stack_fun_LMC_2_downsample(BaseName,Suffix,ImageSize,'PixelCutoff',PixelCutoffRegistration,'XCorrBounds',XCorrBounds,'RegisterColorChannels',1,'NumPar',NumPar,'BeadseqCodePath',BeadseqCodePath,'BarcodeSequence',BarcodeSequence);
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end


%% Bead calling and sequencing
%we now run Bead_Seq itself. Again, where parallelization is trivial in
%Bead_Seq, it is implemented naively using parfor.

%If you have run the basecalling previously, and are rerunning it, you
%should rename that folder. For simplicity, we don't date the basecalling
%folder.

display('Base Calling')


for puck=1:length(PuckNames) %note we are trying to run this overnight without parfor, because we can't figure out why it keeps crashing, and it only crashes for large sized images, and only on puck 8_7
    %We now send the command to do the registration with find_roi
    BaseName=[ProcessedImageFolders{puck},PuckNames{puck},'_Ligation_'];
    suffix='_Stitched';
    disp(['Beginning basecalling for puck number ',num2str(puck)])
	[Bead BeadImage]=BeadSeqFun6_FC_downsample(BaseName,suffix,OutputFolders{puck},BeadZeroThreshold,BarcodeSequence,NumPar,NumLigations,PuckNames{puck},EnforceBaseBalance,BaseBalanceTolerance,'PixelCutoff',PixelCutoffBasecalling,'DropBases',DropBases,'BeadSizeThreshold',BeadSizeCutoff,'PuckImageSubstraction',PuckImageSubstraction);
    
    %added by Jilong for new barcode matching
    [UniqueBeadBarcodes,BBFirstRef,BBOccCounts]=unique([Bead.Barcodes]);
    UniqueBeadLocations={Bead(BBFirstRef).Locations};
    BaseBalanceBarcodes=[UniqueBeadBarcodes];
    BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,NumBases)},'UniformOutput',false);
    BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};
    UniqueBeadBarcodesForExport=char(replace(BaseBalanceBase5Barcodes,{'0','1','2','3','4'},{'N','T','G','C','A'}));
    if NumBases<14 %This is to deal with the InverseLigationSequence -- the export barcodes have to be 14 bases long
        UniqueBeadBarcodesForExport(:,NumBases+1:14)='N';
        UniqueBeadBarcodesForExport=UniqueBeadBarcodesForExport(:,WhichLigationsAreMissing);
    end
    %UniqueBeadBarcodesForExport=UniqueBeadBarcodesForExport(:,InverseLigationSequence);
    file=fullfile(OutputFolders{puck},'BeadBarcodes.txt');
    dlmwrite(file,[UniqueBeadBarcodesForExport]);
    file=fullfile(OutputFolders{puck},'BeadLocations.txt');
    dlmwrite(file,[UniqueBeadLocations]);
    %added done
%The outputted files are of the form 
%[BaseName,int2str(mm),' channel ',int2str(k),suffix,' transform.tif']
end


%% Delete intermediate files

if lower(DeleteIntermediateFiles)~="false"
	for puck=1:length(PuckNames)
		rmdir(ProcessedImageFolders{puck}, 's');
	end
end


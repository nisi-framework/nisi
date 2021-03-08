function zipResultBuiltFiles(nameDataset)
% Dieses Programm schreibt das aktuelle Analysefile und die Kopien der
% Configfiles aller Teildatensätze im nameDataset Ordner in ein Zip-file 
%
% Ursprungsversion von
% Fabian Hufgard
%
% INPUT
%     nameDataset - Welcher Datensatz (im Data-Folder) soll gesichert werden
%
%
% OUTPUT
%     Zip file - Inhalt:
%        Analyze File
%        Config_Copy Files
%
%
% Zukünftig auch Rohdaten dazuspeichern?

% Define directory of nameDataset folder and change current folder to it
load pfade.mat
dirNameDataset = ([pathDataFolder nameDataset]);
cd(dirNameDataset) 

% List folder content 
contentNameDataset = dir;

% Define name of destination zip file
theDate=num2str(yyyymmdd(datetime));
theDate=theDate(3:end);
nameZipFile = ([dirNameDataset '\' theDate '_' nameDataset '_BuildFiles_0.zip']); 
   % Make sure not to overwrite older zip files from the same day
   while exist(nameZipFile,'file') == 2         % as long as destination file name exists, rename file
         n = str2num(nameZipFile(end-4));       % get count of currently checked file
         nameZipFile(end-4) = int2str(n + 1);   % rename file to next number
   end
   
% Create temporary folder
folderTemp = nameZipFile(1:end-4);
mkdir(folderTemp);


% Find analyze file and create list of all files (in case there are more than 1)
for i = 3:length(contentNameDataset)   % starts at 3, because entries 1 and 2 are . and ..
   if ~contentNameDataset(i).isdir && contentNameDataset(i).name(1) == 'a' % If filename starts with a
      nameAF = contentNameDataset(i).name;
      copyfile(nameAF,folderTemp);  %Copy found analyze files to temporary folder
   end
end

% Find subfolders
nSF=0;
for i = 3:length(contentNameDataset) % starts at 3, because entries 1 and 2 are . and ..
   if contentNameDataset(i).isdir
      currFolderFull = [contentNameDataset(i).folder '\' contentNameDataset(i).name];
      if ~strcmp(currFolderFull, folderTemp)
         nSF=nSF+1;
         listSFs{nSF}=contentNameDataset(i).name;
      end
   end
end

% Find config files in subfolders and write to config list
for j = 1:length(listSFs)
   cd([dirNameDataset '/' listSFs{j}]);
   contentCurrentSF = dir; 
   for k = 1:length(contentCurrentSF) % check all files in SF for config copies
      curFN = contentCurrentSF(k).name;
      if ~contentCurrentSF(k).isdir
         if strcmp(curFN(end-13:end),'_Config_Copy.m')
            copyfile(curFN,folderTemp);
         end
      end
   end
end

% Create zip file
cd(folderTemp)
T1 = dir;
for p = 3:length(T1)
   T2{p-2} = T1(p).name; 
end
zip(nameZipFile,T2);

% Delete temporary folder
cd(dirNameDataset);
rmdir(folderTemp,'s');

end

function varargout = Analyser(varargin)
% ANALYSER MATLAB code for Analyser.fig
%      ANALYSER, by itself, creates a new ANALYSER or raises the existing
%      singleton*.
%
%      H = ANALYSER returns the handle to a new ANALYSER or the handle to
%      the existing singleton*.
%
%      ANALYSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSER.M with the given input arguments.
%
%      ANALYSER('Property','Value',...) creates a new ANALYSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Analyser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Analyser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Analyser

% Last Modified by GUIDE v2.5 05-Aug-2017 22:09:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Analyser_OpeningFcn, ...
                   'gui_OutputFcn',  @Analyser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Analyser is made visible.
function Analyser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Analyser (see VARARGIN)
% Flag for clicking the button.

% Choose default command line output for Analyser
handles.output = hObject;

% Condition necessary to have always an option for running mode checked.
if get(handles.offlineModeChkBox,'Value') == 1
    set(handles.offlineModeChkBox,'Enable', 'inactive');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Analyser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Analyser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadFileBtn.
function loadFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% LOAD FILE BUTTON.

% When loading a new file is necessary to "clean" the plots and text views for end user.

% Clean graph of raw data.
cla(handles.figureRawData,'reset');
set(handles.figureRawData,'XMinorTick', 'on');
set(handles.figureRawData,'XGrid', 'on');
set(handles.figureRawData,'YGrid', 'on');
set(handles.figureRawData, 'UIContextMenu', handles.none);

% Clean graph of filtered data.
cla(handles.figureFilters,'reset');
set(handles.figureFilters,'XMinorTick', 'on');
set(handles.figureFilters,'XGrid', 'on');
set(handles.figureFilters,'YGrid', 'on');
set(handles.figureFilters, 'UIContextMenu', handles.none);

% Clean graph of stress zone data.
cla(handles.figureStressDetect,'reset');
set(handles.figureStressDetect,'XMinorTick', 'on');
set(handles.figureStressDetect,'XGrid', 'on');
set(handles.figureStressDetect,'YGrid', 'on');
set(handles.figureStressDetect, 'UIContextMenu', handles.none);

% Reseting values for end user.
set(handles.filterLowPassValue, 'String', '0.00 Hz');
set(handles.filterBandPassValue, 'String', '[ 0.00 - 0.00 ] Hz');
set(handles.stressPeakValue, 'String', '0.00 s');
set(handles.stressStartValue, 'String', '0.00 s');
set(handles.stopStressValue, 'String', '0.00 s');
set(handles.fileName, 'String', 'File name');
set(handles.warningAnalyse, 'Visible', 'on', 'String', 'Load a valid file for analysing, please.');
set(handles.stressAppliedValue, 'String', 'None');
set(handles.uploadDataBtn,'Enable', 'off');

% Managing warnings text views.
set(handles.loadFileInfoValue, 'String', 'Load a valid file, please.');
set(handles.exportDataBtn, 'Enable', 'off');

% Below is the two main option off this software.
% If the offline mode is checked, the user can upload a file from this
% computer normally.
% Else the end user can research a file directly on the database for
% analysing.

%% OFFLINE MODE
if get(handles.offlineModeChkBox,'Value') == 1
    
    % Looking for the exeperiment file
    [FileName, PathName] = uigetfile('*.csv','Select the file');
    if isequal(FileName,0) || isequal(PathName,0)
        set(handles.analyseBtn, 'Enable', 'off');
        set(handles.warningAnalyse, 'Visible', 'on', 'String', 'Load a valid file for analysing, please.');
        return
    end

    % Wait bar
    waitBarLoadFile = waitbar(0, 'Please wait...');
    
    % Larger view setting.
    set(handles.figureRawData, 'UIContextMenu', handles.plotLarger_RawData);

    %Info for saving the data
    [pathStr,nameStr,extStr] = fileparts([PathName FileName]);
    handles.pathStr_H = pathStr;
    handles.nameStr_H = nameStr;
    handles.extStr_H = extStr;
    handles.isUploadAvailable_H = false;

    set(handles.fileName, 'String', nameStr);

    nomDuFichierCSV = [PathName FileName];

    % Looking only for the right part of the file
    fichierCSV = csvread(nomDuFichierCSV, 1, 3);
    % Reading the timespamps of the experiment
    fid = fopen(nomDuFichierCSV);
    dateTime = textscan(fid, '%*s %s %*[^\n]','HeaderLines',1, 'Delimiter', ',');
    fclose(fid);
    % Casting for MATLAB variables.
    dateTime = cellstr(dateTime{1,:});
    heureDebut = dateTime{1,1};
    dateTime = datetime(dateTime);
    duration = seconds(dateTime(length(dateTime)) - dateTime(1));

    nombredeMesures = length(fichierCSV(:,1));

    pause(0.1);
    waitbar(50/100);

    rawData = fichierCSV(:,7:end);

    dataBrute = [];
    for i = 1:nombredeMesures
        dataBrute = horzcat(dataBrute, rawData(i,:));
    end
    
    handles.dataBrute_H = dataBrute;

    fs = length(dataBrute)/duration;
    handles.fs_H = fs;

    t = 1:(duration-1)/length(dataBrute):duration;
    t = t(1:end-(length(t) - length(dataBrute)));
    handles.t_H = t;

    % Wait bar
    waitbar(100/100);
    pause(0.1);
    close(waitBarLoadFile);

else
%% ONLINE MODE
    
    % add the downloaded JAR file to java path
    javaaddpath('mongo-java-driver-2.13.1.jar');

    % import the library
    import com.mongodb.*;
    
    db = handles.db_H;
    m = handles.m_H;
    
    % Wait bar
    waitBarLoadFile = waitbar(0, 'Please wait... Getting data from database.');
    
    % Getting the name of the file and verifying if the end user texted
    % something in the edit test.
    if get(handles.tapFileNameDataBaseValue, 'Value') == 1
        nameStr = get(handles.tapFileNameDataBaseValue, 'String');
    else
        pause(0.1);
        close(waitBarLoadFile);
        return
    end

    % get a collection
    ensemble_donnes_collection = db.getCollection('ensemble_donnees');

    % finding the experience
    dbObjTestEnsembleDonnees = BasicDBObject('nom_fichier_csv', nameStr);
    myDoc_ensemble_donnes_collection = ensemble_donnes_collection.findOne(dbObjTestEnsembleDonnees);
    
    % Verifying id the raw data is already online
    if  isempty(myDoc_ensemble_donnes_collection)
        pause(0.1);
        close(waitBarLoadFile);
        return
    end
    
    % get the the number of data per line
    dataPerLigne = myDoc_ensemble_donnes_collection.get('frequence_utilisee');
    
    % Copy of the DB Object for getting other features necessary for
    % uploading the analysed data correctly.
    dbObjTestCopyEnsembleDonnees = ensemble_donnes_collection.findOne(dbObjTestEnsembleDonnees).copy();

    % get the indexof the objects for using in the other collection.
    releves_donnees = myDoc_ensemble_donnes_collection.get('id_releves_donnees').toArray;
    releves_donnees = cell2mat(cell(releves_donnees))';

    % get a collection
    releve_donnees_collection = db.getCollection('releve_donnees');

    % Verifying if the data is already analysed and uploaded.
    % get a collection
    filtre_passe_bas_collection = db.getCollection('aaa_filtre_passe_bas');
    % get a collection
    moyennage_collection = db.getCollection('aab_moyennage');

    waitbar(15/100);

    %% VERIFYING IF THE ANALYSED DATA IS ALREADY ONLINE
    
    dbObjTestIfExist = BasicDBObject('id_', releves_donnees(1));
    iddbObjTestIfExist = filtre_passe_bas_collection.findOne(dbObjTestIfExist);
    if  isempty(iddbObjTestIfExist)
        handles.isUploadAvailable_H = true;
    else
        handles.isUploadAvailable_H = false;
    end

    %% GETTING DATA FROM DATABASE
    
    % allocation of vectors
    dataBrute = []; % vector of raw data.
    liste_date_complet = []; % vector for all timestamps of the experience. 

    % as the data was split for uploading from the .csv we need to "rebuild"
    % the data in order to analyse.
    qqte_donnees_par_objet = []; % vector for how many lines of data has each object.

    % get the all the raw data and put it in Matlab variables.
    for i = 1 : length(releves_donnees)
        dbObjTestReleveDonnees = BasicDBObject('id_', releves_donnees(i));
        myDoc_releve_donnees_collection = releve_donnees_collection.findOne(dbObjTestReleveDonnees).keySet.toArray;
        liste_dataTime = cell(myDoc_releve_donnees_collection);
        liste_dataTime_char = char(liste_dataTime(4:end, :));
        liste_date_complet = vertcat(liste_date_complet, liste_dataTime_char);

        objet = releve_donnees_collection.findOne(dbObjTestReleveDonnees);
        qqte_donnees_par_objet = [qqte_donnees_par_objet length(liste_dataTime_char)];

        for j = 1 : length(liste_dataTime_char)
            donnesBrute_par_ligne = objet.get(liste_dataTime_char(j,:)).toArray;
            donnesBrute_par_ligne = str2double(cell(donnesBrute_par_ligne))';

            %Recupere toutes les donnees en un vecteur.
            dataBrute = horzcat(dataBrute, donnesBrute_par_ligne);
        end

        % get the time of the start and end of the experience.
        if i == 1
            heure_debut = datetime(liste_dataTime(4, 1));
        end

        if i == length(releves_donnees)
            heure_fin = datetime(liste_dataTime(length(liste_dataTime), 1));
        end
        waitbar((50+i)/100);
    end
    
    % Wait bar
    waitbar(100/100);
    pause(0.1);
    close(waitBarLoadFile);
    
    duration = seconds(heure_fin - heure_debut);
    fs = length(dataBrute)/duration;
  
    t = 1:(duration-1)/length(dataBrute):duration;
    t = t(1:end-(length(t) - length(dataBrute)));
    
    % Handles
    handles.qqte_donnees_par_objet_H = qqte_donnees_par_objet;
    handles.releves_donnees_H = releves_donnees;
    handles.dbObjTestCopyEnsembleDonnees_H = dbObjTestCopyEnsembleDonnees;
    handles.liste_date_complet_H = liste_date_complet;
    handles.filtre_passe_bas_collection_H = filtre_passe_bas_collection;
    handles.moyennage_collection_H = moyennage_collection;
    handles.dataPerLigne_H = dataPerLigne;
    
    handles.nameStr_H = nameStr(1:end-4);
    handles.t_H = t;
    handles.fs_H = fs;
    handles.dataBrute_H = dataBrute;
    
end

% Flag for clicking the button.
set(handles.warningAnalyse, 'String', 'Analysis can be done ...');
set(handles.analyseBtn, 'Enable', 'on');

% Change text viem for user.
set(handles.loadFileInfoValue, 'String', 'File loaded sucessfully.');

% Update handles structure
guidata(hObject, handles);

% Plot of the raw data.
figureRawData = handles.figureRawData;
axes(handles.figureRawData)
plot(figureRawData, t, dataBrute, 'b');
set(figureRawData,'Xlabel', xlabel('Time (s)'));
set(figureRawData,'Ylabel', ylabel('Amplitude (V)'));
set(figureRawData,'XMinorTick', 'on');
set(figureRawData,'XGrid', 'on');
set(figureRawData,'YGrid', 'on');

%% END OF THE LOAD FILE BUTTON.


% --- Executes on button press in analyseBtn.
function analyseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to analyseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Analysing button.

% Wait bar
waitBarLoadFile = waitbar(0, 'Please wait...');

% Getting handles.
dataBrute = handles.dataBrute_H;
fs = handles.fs_H;
t = handles.t_H;

%% Spectral analysis
Order = 8;

dataBruteFFT = fft(dataBrute);
dataBruteFFTConj = dataBruteFFT.*conj(dataBruteFFT)/length(dataBruteFFT);
dataBruteDB = 20.*log10(dataBruteFFTConj); % Its an amplitude so : AdB = 20.*log10(x)
f = fs/length(dataBruteFFT)*(0:(length(dataBruteFFT)-1));

% Semi global variables
handles.dataBrute_HFFTConj_H = dataBruteFFTConj;
handles.dataBruteDB_H = dataBruteDB;
handles.f_H = f;

%%  Finding la bonne frequence à analyser.
% From a set value ahead the software will look for the max amplitude in the
% spectral analyssis. In this case it is set for 3 Hz.
StarAnalysisFrequency = 3*round(1/(fs/length(dataBruteFFT)));

% Gettin the max amplitude.
[maxAmpl, indexMaxAmplFreq] = max(dataBruteDB(StarAnalysisFrequency:round(length(dataBruteDB)/2)));
indexMaxAmplFreq = indexMaxAmplFreq + StarAnalysisFrequency;

% Vector of 1Hz after the max aplitude
highLim = dataBruteDB(indexMaxAmplFreq:indexMaxAmplFreq+round(1/(fs/length(dataBruteFFT))));

% Getting the high limit of the filter band pass.
dataBruteDBMeanH = smooth(highLim,round(1/(fs/length(dataBruteFFT))*0.2))'; % 20%
indexMinAmplLimHigh = length(highLim(dataBruteDBMeanH > -15)); % 15dB = -15dB = 0.17
if indexMinAmplLimHigh == length(highLim)
    [minAmplLimHigh, indexMinAmplLimHigh] = min(dataBruteDBMeanH);
end

% Vector of 1Hz before the max aplitude
lowLim = dataBruteDB(indexMaxAmplFreq-round(1/(fs/length(dataBruteFFT))):indexMaxAmplFreq);
dataBruteDBMeanL = smooth(lowLim,round(1/(fs/length(dataBruteFFT))*0.2))';% 20%
indexMinAmplLimLow = length(lowLim(dataBruteDBMeanL > -15)); % -15dB = 0.17
if indexMinAmplLimLow == length(lowLim)
    [minAmplLimLow, indexMinAmplLimLow] = min(dataBruteDBMeanL);
end

% Frequency of max amplitude
fBandPass = f(indexMaxAmplFreq);

% Getting th cut frequency of the filters.
fBandPassLimHigh = f(indexMaxAmplFreq + indexMinAmplLimHigh);
fBandPassLimLow = f(indexMaxAmplFreq - (length(dataBruteDBMeanL) - indexMinAmplLimLow));

waitbar(25/100);
pause(0.1);

%% Filtering

bpFilt = designfilt('bandpassiir','FilterOrder',Order, ...
    'HalfPowerFrequency1',fBandPassLimLow, ...
    'HalfPowerFrequency2',fBandPassLimHigh, ...
    'SampleRate',fs);

dataBandFiltre = filter(bpFilt, dataBrute);

lpFilt = designfilt('lowpassiir','FilterOrder',Order, ...
    'PassbandFrequency',fBandPassLimHigh,'PassbandRipple',0.2, ...
    'SampleRate',fs);

 dataLowFiltre = filter(lpFilt, dataBrute);
 
 waitbar(50/100);
 pause(0.1);
                                        
%% STRESS ANALYSIS - IDENTIFICATION

% Getting only the postive part of the filtered band-pass signal.
positifDataFiltre = dataBandFiltre;
for i = 1 : length(dataBandFiltre)
    if dataBandFiltre(i) < 0
    positifDataFiltre(i) = 0;
    end
end

% Vector of one second value according to the signal sample frequency.
unSecond = round(1/((t(end)-t(1))/length(t)));
handles.unSecond_H = unSecond;

% Getting the envelope of the filtered band-pass signal using an
% envelope each 0.5 s.
% In the begin of the signal. There is a issue with the damping of
% the filter. In this case is set for 4s.
dampingtTime = 4;
[envHigh, envLow] = envelope(positifDataFiltre(t > dampingtTime),round(0.5*unSecond),'peak');
[maxAmplTest, indexMaxAmplTest] = max(envHigh(1:end-dampingtTime*unSecond));

% Getting the max of the envelope of the filtered data. 
% Wich is probably the moment of stress in the plant.
[envHighTempDebut, envLowTempDebut] = envelope(positifDataFiltre(t <= dampingtTime),round(0.5*unSecond),'peak');

% Temporary vector for compensate the time of damping.
% Only for fullfil the begin of the main analysis vector.
indexAuxTempDebut = length(t(t <= dampingtTime));
regionOfStressTempDebut = zeros(1, indexAuxTempDebut);

% Getting the stress peak index and relating with time vector.
indexStressPeak = indexMaxAmplTest;
tStressPeak = t(indexMaxAmplTest + indexAuxTempDebut);

% Normalising the difference between the max and the mean of the
% envelope.
difference = maxAmplTest - mean(envHigh);

%% Mean both sides of the peak for getting the start and the begin of the zone stress.

moyenneApres = sum(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond))/length(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond));
moyenneAvant = sum(envHigh(indexMaxAmplTest-unSecond:indexMaxAmplTest))/length(envHigh(indexMaxAmplTest-unSecond:indexMaxAmplTest));

% Getting a parameter normalized for each signal.
differenceN = moyenneAvant - mean(envHigh);

limLow = unSecond;
while differenceN > 0.4*difference
    moyenneAvant = sum(envHigh(indexMaxAmplTest-unSecond-limLow:indexMaxAmplTest))/length(envHigh(indexMaxAmplTest-unSecond-limLow:indexMaxAmplTest));
    differenceN = moyenneAvant - mean(envHigh);
    limLow = limLow + unSecond;
end

% Releasing the difference for the other side.
differenceN = moyenneApres - mean(envHigh);

limHigh = unSecond;
while differenceN > 0.4*difference
    moyenneApres = sum(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond+limHigh))/length(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond+limHigh));
    differenceN = moyenneApres - mean(envHigh);
    limHigh = limHigh + unSecond;
end

% Moment of start and end of the zone stress.
tDebutRegionStress = t(indexMaxAmplTest + indexAuxTempDebut - limLow);
tFinRegionStress = t(indexMaxAmplTest + indexAuxTempDebut + limHigh);

% Designing the stress zone vector.
regionOfStress = ones(1, length(envHigh));

for i = 1 : length(regionOfStress)
    if indexStressPeak - limLow > i
        regionOfStress(i) = 0;
    end
    if indexStressPeak + limHigh < i
        regionOfStress(i) = 0;
    end
end

regionOfStressNormalise = max(envHigh).*regionOfStress;

waitbar(75/100);
pause(0.1);

%% Code block for testing for identification and characterisation of the signal.

% Getting the diff of the envelope.
dataPeak = [0 diff(envHigh)];

dataPeakTemp = dataPeak.*regionOfStress;
dataPeakTemp = [regionOfStressTempDebut dataPeakTemp];
dataPeak = [regionOfStressTempDebut dataPeak];

% Detecting the borders of the stress peak. We get the min and max diff
% of the stress zone, so it means the "real" start and end of the
% stress.
[maxAmplPeak, indexMaxAmplPeak] = max(dataPeakTemp);
[minAmplPeak, indexMinAmplPeak] = min(dataPeakTemp);
firstPeakMax = t(indexMaxAmplPeak);
firstPeakMin = t(indexMinAmplPeak);

% We need the two biggest diffs of the stress zone. So we detect the
% biggest and then we force it to zero in order to detect the next
% biggest diff.
contAuxMaxAmplPeak = indexMaxAmplPeak;
auxSimetric = 0;
while dataPeak(contAuxMaxAmplPeak + auxSimetric) >= 0
    dataPeakTemp(contAuxMaxAmplPeak + auxSimetric) = 0;
    dataPeakTemp(contAuxMaxAmplPeak - auxSimetric) = 0;
    auxSimetric = auxSimetric + 1;
end

contAuxMinAmplPeak = indexMinAmplPeak;
auxSimetric = 0;
while dataPeak(contAuxMinAmplPeak + auxSimetric) <= 0
    dataPeakTemp(contAuxMinAmplPeak + auxSimetric) = 0;
    dataPeakTemp(contAuxMinAmplPeak - auxSimetric) = 0;
    auxSimetric = auxSimetric + 1;
end

[secondMaxAmplPeak, indexSecondMaxAmplPeak] = max(dataPeakTemp);
[secondMinAmplPeak, indexSecondMinAmplPeak] = min(dataPeakTemp);
segundoPicoMax = t(indexSecondMaxAmplPeak);
segundoPicoMin = t(indexSecondMinAmplPeak);

% If the second max is lower then 50% of the first one the second peak is
% changed for value of the first one.
if abs(secondMaxAmplPeak) < 0.5*abs(maxAmplPeak)
    secondMaxAmplPeak = maxAmplPeak;
    indexSecondMaxAmplPeak = indexMaxAmplPeak;
end

if abs(secondMinAmplPeak) < 0.5*abs(minAmplPeak)
    secondMinAmplPeak = minAmplPeak;
    indexSecondMinAmplPeak = indexMinAmplPeak;
end

amplitudesVector = [maxAmplPeak minAmplPeak secondMaxAmplPeak secondMinAmplPeak];
indexAmplitudesVector = [indexMaxAmplPeak indexMinAmplPeak indexSecondMaxAmplPeak indexSecondMinAmplPeak];
indexAmplitudesVector = indexAmplitudesVector(amplitudesVector ~= 0);

% Vector temporary only for showing a nice view for end user.
newPeriodStress = zeros(1, length(dataPeak));

maxDataPeak = max(dataPeak);
maxEnvHigh = max(envHigh);
minIndexAmplitudesVector = min(indexAmplitudesVector);
maxIndexAmplitudesVector = max(indexAmplitudesVector);

for j = 1 : length(newPeriodStress)
    if minIndexAmplitudesVector < j && maxIndexAmplitudesVector > j
        newPeriodStress(j) = maxEnvHigh;
    end
end

% Relating the diff vector with the "real" stress vector. In order to
% determine where exctly is the stress.
limNeg = 0;
limPos = 0;
for w = 2: length(dataPeak)
    if newPeriodStress(w) == maxDataPeak
        for n = 1 : w-1
            if limPos == 1 && limNeg == 1
                break
            end
            if w-n > 0
                if abs(dataPeak(w-n)) > 0.01*maxDataPeak && limNeg == 0
                    newPeriodStress(w) = maxDataPeak;
                    newPeriodStress(w-n) = maxDataPeak;
                else
                    limNeg = 1;
                end
            end
            if w+n < length(newPeriodStress)
                if abs(dataPeak(w+n)) > 0.01*maxDataPeak && limPos == 0
                    newPeriodStress(w) = maxDataPeak;
                    newPeriodStress(w+n) = maxDataPeak;
                else
                    limPos = 1;
                end
            end
        end
        limNeg = 0;
        limPos = 0;
    end
end

% Getting the "real" start and end of the stress. 
isInverted = 0;
for p = 2 : length(newPeriodStress)
    if newPeriodStress(p) ~= newPeriodStress(p-1)
        if isInverted == 0
            tStressStart = t(p);
            isInverted = isInverted + 1;
        end
        if isInverted == 1
            tStressEnd = t(p);
        end
    end
end

if length(newPeriodStress(newPeriodStress ~= 0)) > 5*unSecond
    set(handles.stressAppliedValue, 'String', 'Flame')
else
    set(handles.stressAppliedValue, 'String', 'Cut or touch')
end

%% ENABLING/DISABLING UPLOAD BUTTON IF THE ANALYSED DATA IS ALREADY ONLINE

if  handles.isUploadAvailable_H
    set(handles.uploadDataBtn,'Enable', 'on');
else
    set(handles.uploadDataBtn,'Enable', 'off');
end

%% Semi global variables

handles.dataLowFiltre_H = dataLowFiltre;
handles.dataBandFiltre_H = dataBandFiltre;
handles.regionOfStressNormalise_H = [regionOfStressTempDebut regionOfStressNormalise];
handles.envelopeFilteredSignal_H = [envHighTempDebut envHigh];
handles.fBandPassLimHigh_H = fBandPassLimHigh;
handles.fBandPassLimLow_H = fBandPassLimLow;
handles.fBandPass_H = fBandPass;
handles.lowLim_H = lowLim;
handles.highLim_H = highLim;
handles.indexMaxAmplFreq_H = indexMaxAmplFreq;
handles.dataBruteDBMeanL_H = dataBruteDBMeanL;
handles.dataBruteDBMeanH_H = dataBruteDBMeanH;

%% Setting values for end user.
set(handles.filterLowPassValue, 'String', [num2str(fBandPassLimHigh) ' Hz']);
set(handles.filterBandPassValue, 'String', ['[ ' num2str(fBandPassLimLow) '  -  ', ...
                                            num2str(fBandPassLimHigh) ' ] Hz']);
set(handles.stressPeakValue, 'String', [num2str(tStressPeak) ' s']);
set(handles.stressStartValue, 'String', [num2str(tStressStart) ' s']);
set(handles.stopStressValue, 'String', [num2str(tStressEnd) ' s']);
set(handles.figureFilters, 'UIContextMenu', handles.plotLarger_DataFiltered);
set(handles.figureStressDetect, 'UIContextMenu', handles.plotLarger_StressDetect);

waitbar(100/100);
pause(0.1);
close(waitBarLoadFile);

%% Plotting figures of the filter reponses.

% Plot of the filtered data.
figureFilters = handles.figureFilters;
axes(handles.figureFilters)

% Low-pass
plot(figureFilters, t, dataLowFiltre, 'g');
set(figureFilters,'Xlabel', xlabel('Time (s)'));
set(figureFilters,'Ylabel', ylabel('Amplitude (V)'));
set(figureFilters,'XMinorTick', 'on');
set(figureFilters,'XGrid', 'on');
set(figureFilters,'YGrid', 'on');
hold on

% Low-pass
plot(figureFilters, t, dataBandFiltre, 'b');
hold off

% Stress analysis
figureStressDetect = handles.figureStressDetect;
axes(handles.figureStressDetect)

plot(figureStressDetect, t, [envHighTempDebut envHigh], 'b');
set(figureStressDetect,'Xlabel', xlabel('Time (s)'));
set(figureStressDetect,'Ylabel', ylabel('Amplitude (V)'));
set(figureStressDetect,'XMinorTick', 'on');
set(figureStressDetect,'XGrid', 'on');
set(figureStressDetect,'YGrid', 'on');
hold on

plot(figureStressDetect, t, [regionOfStressTempDebut regionOfStressNormalise], 'r');
hold off

%% Managing warnings.
set(handles.warningAnalyse, 'String', 'Analysis DONE!');
pause(1);
set(handles.analyseBtn, 'Enable', 'off');
set(handles.warningAnalyse, 'String', 'The data can be exported.');
set(handles.exportDataBtn, 'Enable', 'on');

%% Update handles structure
guidata(hObject, handles);
%% END OF THE BUTTON ANALYSE.


% --- Executes on button press in exportDataBtn.
function exportDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to exportDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

directoryName= uigetdir('C:\', 'Select a directory for saving the data');
if directoryName == 0
    set(handles.warningAnalyse, 'String', 'Select a folder for exporting data, please.');
    return
end

% Getting the name of the file
nameStr = handles.nameStr_H;

figure('units','normalized','outerposition',[0 0 1 1])

%% Raw data

% Semi global variables
dataBrute = handles.dataBrute_H;
t = handles.t_H;

% Saving the data
save([directoryName '\' nameStr '_time_vector'], 't');
save([directoryName '\' nameStr '_raw_data_vector'],'dataBrute');

% Save plot
plot(t, dataBrute, 'b');
grid on
title(['Raw data for ' nameStr])
xlabel('Time (s)')
ylabel('Amplitude (V)')
legend('Raw data')
saveas(gcf, fullfile(directoryName, [nameStr '_raw_data']), 'jpeg');

%% Filtered data

% Semi global variables
dataLowFiltre = handles.dataLowFiltre_H; 
dataBandFiltre = handles.dataBandFiltre_H;
regionOfStressNormalise = handles.regionOfStressNormalise_H;
envelopeFilteredSignal = handles.envelopeFilteredSignal_H;
fBandPassLimHigh = handles.fBandPassLimHigh_H;
fBandPassLimLow = handles.fBandPassLimLow_H;

% Save plot
plot(t, dataLowFiltre, 'g');
grid on
title(['Low-pass filtered signal for ' nameStr])
xlabel('Time (s)')
ylabel('Amplitude (V)')
legend(['Low-pass filter signal at ' num2str(fBandPassLimHigh) ' Hz'])
saveas(gcf, fullfile(directoryName, [nameStr '_low_pass_filtered_signal']), 'jpeg');

plot(t, dataBandFiltre, 'b');
grid on
title(['Band-pass filtered signal for ' nameStr])
xlabel('Time (s)')
ylabel('Amplitude (V)')
legend(['Band-pass filter between [ ' num2str(fBandPassLimLow) ' - ' num2str(fBandPassLimHigh) ' ] Hz'])
saveas(gcf, fullfile(directoryName, [nameStr '_band_pass_filtered_signal']), 'jpeg');

% Saving the data
save([directoryName '\' nameStr '_low_pass_filtered_signal_vector'], 'dataLowFiltre');
save([directoryName '\' nameStr '_band_pass_filtered_signal_vector'],'dataBandFiltre');

%% Analyse Spectrale

% Semi global variables
dataBruteFFTConj = handles.dataBrute_HFFTConj_H;
dataBruteDB = handles.dataBruteDB_H;
f = handles.f_H;
fBandPass = handles.fBandPass_H;
lowLim = handles.lowLim_H;
highLim = handles.highLim_H;
indexMaxAmplFreq = handles.indexMaxAmplFreq_H;
dataBruteDBMeanL = handles.dataBruteDBMeanL_H;
dataBruteDBMeanH = handles.dataBruteDBMeanH_H;

% Saving the data
save([directoryName '\' nameStr '_spectral_analysys_vector'], 'dataBruteFFTConj');
save([directoryName '\' nameStr '_spectral_analysys_dB_vector'],'dataBruteDB');

%Saving and Plotting detail of the frequency filtered
plot(f(indexMaxAmplFreq - length(lowLim):indexMaxAmplFreq + length(highLim) - 1), [lowLim highLim])
hold on
plot(f(indexMaxAmplFreq - length(lowLim):indexMaxAmplFreq + length(highLim) - 1), [dataBruteDBMeanL dataBruteDBMeanH])
grid on
title(['Spectral analysis in detail for ' nameStr])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend(['Frequency of max amplitude is ' num2str(fBandPass) ' Hz'], 'Smooth');
saveas(gcf, fullfile(directoryName, [nameStr '_spectral_analysis_detail']), 'jpeg');

% Save plot
subplot(2, 1, 1)
plot(f(1:(round(length(f)/10))), dataBruteDB(1:(round(length(dataBruteDB)/10))), 'b');
title(['Spectral analysis of raw data for ' nameStr])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
subplot(2, 1, 2)
plot(f(1:(round(length(f)/10))), dataBruteFFTConj(1:(round(length(dataBruteFFTConj)/10))), 'b');
xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on
saveas(gcf, fullfile(directoryName, [nameStr '_spectral_analysis']), 'jpeg');

%% Closing the figure
close;

%% Finishing export data
set(handles.exportDataBtn, 'Enable', 'off');
set(handles.warningAnalyse, 'String', 'Data exported!');
%%END OF THE BUTTON EXPORT DATA


% --- Displays nothing at larger size in a new figure.
function none_Callback(hObject, eventdata, handles)
% hObject    handle to none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Displays contents of figureRawData at larger size in a new figure.
function plotLarger_RawData_PopUp_Callback(hObject, eventdata, handles)
% hObject    handle to plotLarger_RawData_PopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Getting the name of the file
nameStr = handles.nameStr_H;

% Create a figure to receive this axes' data
figureRawDataExtFig = figure('units','normalized','outerposition',[0 0 1 1]);
% Copy the axes and size it to the figure
figureRawDataCopy = copyobj(handles.figureRawData,figureRawDataExtFig);
set(figureRawDataCopy,'Units','Normalized',...
              'Position',[.05,.20,.90,.60]);
% Assemble a title for this new figure
title(['Raw data for ' nameStr],'Fontweight','bold')
legend('Raw data')


% --- Displays contents of figureFilters at larger size in a new figure.
function plotLarger_DataFiltered_PopUp_Callback(hObject, eventdata, handles)
% hObject    handle to plotLarger_DataFiltered_PopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Getting the settings of the file
nameStr = handles.nameStr_H;
fBandPassLimHigh = handles.fBandPassLimHigh_H;
fBandPassLimLow = handles.fBandPassLimLow_H;

% Create a figure to receive this axes' data
figureFiltersExtFig = figure('units','normalized','outerposition',[0 0 1 1]);
% Copy the axes and size it to the figure
figureFiltersCopy = copyobj(handles.figureFilters, figureFiltersExtFig);
set(figureFiltersCopy,'Units','Normalized',...
              'Position',[.05,.20,.90,.60]);
% Assemble a title for this new figure
title(['Filtered signal for ' nameStr],'Fontweight','bold')
legend(['Low-pass filter at ' num2str(fBandPassLimHigh) ' Hz'] ,...
    ['Band-pass filter between [ ' num2str(fBandPassLimLow) ' - ' num2str(fBandPassLimHigh) ' ] Hz']);


% --- Displays contents of figureStressDetect at larger size in a new figure.
function plotLarger_StressDetected_PopUp_Callback(hObject, eventdata, handles)
% hObject    handle to plotLarger_StressDetected_PopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Getting the name of the file
nameStr = handles.nameStr_H;

% Create a figure to receive this axes' data
figureStressDetectExtFig = figure('units','normalized','outerposition',[0 0 1 1]);
% Copy the axes and size it to the figure
figureStressDetectCopy = copyobj(handles.figureStressDetect,figureStressDetectExtFig);
set(figureStressDetectCopy,'Units','Normalized',...
              'Position',[.05,.20,.90,.60]);
% Assemble a title for this new figure
title(['Region of stress for ' nameStr],'Fontweight','bold')
legend('Filtered signal', 'Stress zone')


function lowPassFilterLegend_Callback(hObject, eventdata, handles)
% hObject    handle to lowPassFilterLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowPassFilterLegend as text
%        str2double(get(hObject,'String')) returns contents of lowPassFilterLegend as a double

% --- Executes during object creation, after setting all properties.
function lowPassFilterLegend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowPassFilterLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bandPassFilterLegend_Callback(hObject, eventdata, handles)
% hObject    handle to bandPassFilterLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bandPassFilterLegend as text
%        str2double(get(hObject,'String')) returns contents of bandPassFilterLegend as a double

% --- Executes during object creation, after setting all properties.
function bandPassFilterLegend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bandPassFilterLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function signalFilteredLegend_Callback(hObject, eventdata, handles)
% hObject    handle to signalFilteredLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of signalFilteredLegend as text
%        str2double(get(hObject,'String')) returns contents of signalFilteredLegend as a double

% --- Executes during object creation, after setting all properties.
function signalFilteredLegend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to signalFilteredLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function stressZoneLegend_Callback(hObject, eventdata, handles)
% hObject    handle to stressZoneLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stressZoneLegend as text
%        str2double(get(hObject,'String')) returns contents of stressZoneLegend as a double

% --- Executes during object creation, after setting all properties.
function stressZoneLegend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stressZoneLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function plotLarger_RawData_Callback(hObject, eventdata, handles)
% hObject    handle to plotLarger_RawData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plotLarger_DataFiltered_Callback(hObject, eventdata, handles)
% hObject    handle to plotLarger_DataFiltered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plotLarger_StressDetect_Callback(hObject, eventdata, handles)
% hObject    handle to plotLarger_StressDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on connect button press for connecting to the database.
function connectDbBtn_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to connectDbBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% INPUT DATABASE INFOS
prompt = {'What is the DATABASE name? (ex. exampleMyDb)',
'What is the IP adress of the database? (ex. 127.0.0.0)',
'What is the PORT for connecting to the database? (ex. 27017)'};
dlg_title = 'Database infos';
sizeOfTheBox = [1 70];
defaultans = {' ', ' ', ' '};
answer = inputdlg(prompt, dlg_title, sizeOfTheBox, defaultans);

if isequal(answer, defaultans) || isempty(answer)
    return
end
dataBaseName = answer{1};
dataBaseIP = answer{2};
dataBasePORT = answer{3};

% Check if all are fullfilled
for i = 1: length(answer)
    if answer{i} == ' '
        return
    end
end

[login, password] = logindlg(); % Getting login and password.
if isempty(login) || isempty(password)
    return
end

%% CONNECTING TO THE DATABASE
% in order to communicate to the DB (mongoDB) we need to use some commands
% from Java language, so we need to put the driver and import the library
% in our path.

% add the downloaded JAR file to java path
javaaddpath('mongo-java-driver-2.13.1.jar');

% import the library
import com.mongodb.*;

% Connect to remote server
m = Mongo(dataBaseIP, str2num(dataBasePORT));

% open DB
db = m.getDB(dataBaseName);

% Authentication
auth = db.authenticate(login , password);
if auth
    set(handles.dataBaseConnectionStatusValue, 'String', 'ONLINE');
else
    set(handles.dataBaseConnectionStatusValue, 'String', 'NOT CONNECTED');
    pause(1);
    set(handles.dataBaseConnectionStatusValue, 'String', 'OFFLINE');
    return
end

handles.m_H = m;
handles.db_H = db;

set(handles.connectDbBtn, 'Enable', 'off');
set(handles.dataBaseConnectionStatusValue, 'Value', 1);
set(handles.disconnectDbBtn, 'Enable', 'on');
set(handles.onlineModeChkBox,'Enable', 'on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on disconnect button press if the user is connect to the database.
function disconnectDbBtn_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to disconnectDbBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% add the downloaded JAR file to java path
javaaddpath('mongo-java-driver-2.13.1.jar');

% import the library
import com.mongodb.*;

m = handles.m_H;

m.close;
m.close();
close(m);

%% Setting infos for end user
set(handles.connectDbBtn, 'Enable', 'on');
set(handles.disconnectDbBtn, 'Enable', 'off');
set(handles.uploadDataBtn, 'Enable', 'off');
set(handles.dataBaseConnectionStatusValue, 'String', 'OFFLINE');
set(handles.dataBaseConnectionStatusValue, 'Value', 0);
set(handles.onlineModeChkBox,'Value', 0);
set(handles.offlineModeChkBox,'Value', 1);
set(handles.onlineModeChkBox,'Enable', 'off');
set(handles.tapFileNameDataBase,'Enable', 'off');
set(handles.tapFileNameDataBaseValue,'Enable', 'off');

% Update handles structure
guidata(hObject, handles);


function tapFileNameDataBaseValue_Callback(hObject, eventdata, handles)
% hObject    handle to tapFileNameDataBaseValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tapFileNameDataBaseValue as text
%        str2double(get(hObject,'String')) returns contents of tapFileNameDataBaseValue as a double
if get(hObject,'Value') == 0
    set(hObject,'Value', 1);
end


% --- Executes during object creation, after setting all properties.
function tapFileNameDataBaseValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tapFileNameDataBaseValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in onlineModeChkBox.
function onlineModeChkBox_Callback(hObject, eventdata, handles)
% hObject    handle to onlineModeChkBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onlineModeChkBox

if get(hObject,'Value') == 1
    set(handles.offlineModeChkBox,'Value', 0);
    set(hObject,'Enable', 'inactive');
    
    set(handles.tapFileNameDataBase,'Enable', 'on');
    set(handles.tapFileNameDataBaseValue,'String', 'Tap the name of the file');
    set(handles.tapFileNameDataBaseValue,'Enable', 'on');
end

set(handles.offlineModeChkBox,'Enable', 'on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in offlineModeChkBox.
function offlineModeChkBox_Callback(hObject, eventdata, handles)
% hObject    handle to offlineModeChkBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of offlineModeChkBox

if get(hObject,'Value') == 1
    set(handles.onlineModeChkBox,'Value', 0);
    set(hObject,'Enable', 'inactive');
    
    set(handles.tapFileNameDataBase,'Enable', 'off');
    set(handles.tapFileNameDataBaseValue,'String', 'Tap the name of the file');
    set(handles.tapFileNameDataBaseValue,'Enable', 'off');
end

set(handles.onlineModeChkBox,'Enable', 'on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on upload button press if the user is connect to the database.
function uploadDataBtn_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uploadDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Wait bar
waitBarUploadFile = waitbar(0, 'Please wait... Uploading data...');

% add the downloaded JAR file to java path
javaaddpath('mongo-java-driver-2.13.1.jar');

% import the library
import com.mongodb.*;

% Getting handles...
qqte_donnees_par_objet = handles.qqte_donnees_par_objet_H;
releves_donnees = handles.releves_donnees_H;
dbObjTestCopyEnsembleDonnees = handles.dbObjTestCopyEnsembleDonnees_H;
liste_date_complet = handles.liste_date_complet_H;
filtre_passe_bas_collection = handles.filtre_passe_bas_collection_H;
moyennage_collection = handles.moyennage_collection_H;
dataPerLigne = handles.dataPerLigne_H;

dataLowFiltre = handles.dataLowFiltre_H; 
dataBandFiltre = handles.dataBandFiltre_H;


%% Upload de donnees vers la BDD.

% fist package of the signal processing - FILTRE PASSE BAS.

% vector of DBObjects.
innerDoc1_Objects = {};

% there is a limite concerning the size of the file upload, it is :
% 16793600 bit. So just like for uploading the file .csv in the DB we need to
% split the data too. It's is usefull because it keeps the same structure
% of other collections too (see releve_donnees).

num_experience = int32(dbObjTestCopyEnsembleDonnees.get('id_'));
data = 1;
k = 0;
for i = 1:length(qqte_donnees_par_objet)
    innerDoc1 = BasicDBObject();
    innerDoc1.put('id_', int32(releves_donnees(i)));
    innerDoc1.put('id_ensemble_donnees', num_experience);
    for j = 1 : qqte_donnees_par_objet(i)
        innerDoc1.put(liste_date_complet(k + j,:), dataLowFiltre(data:(data + dataPerLigne - 1)));
        data = data + dataPerLigne -1;
    end
    k = k + qqte_donnees_par_objet(i);
    innerDoc1_Objects = [innerDoc1_Objects innerDoc1];
end

waitbar(50/100);

% second package of the signal processing - MOYENNAGE.

% vector of DBObjects.
innerDoc2_Objects = {};

data = 1;
k = 0;
for i = 1:length(qqte_donnees_par_objet)
    innerDoc2 = BasicDBObject();
    innerDoc2.put('id_', int32(releves_donnees(i)));
    innerDoc2.put('id_ensemble_donnees', num_experience);
    for j = 1 : qqte_donnees_par_objet(i)
        innerDoc2.put(liste_date_complet(k + j,:), dataBandFiltre(data:(data + dataPerLigne - 1)));
        data = data + dataPerLigne -1;
    end
    k = k + qqte_donnees_par_objet(i);
    innerDoc2_Objects = [innerDoc2_Objects innerDoc2];
end

% upload the data to the DB using the command insert().
for j = 1 : length(qqte_donnees_par_objet)
    
    % Wait for acknowledgement, but don't wait for secondaries to replicate == SAFE
    % insert the document into the database. If the database does not exist
    % and/or if the collection does not exist, it is created at the first
    % insertion. See mongoDB documentation for details.
    wc1 = com.mongodb.WriteConcern(1);
    wc2 = com.mongodb.WriteConcern(1);
    
    % The insert method takes a collection of documents as input, so we must
    % tell the server the writeconcern to control how to write in the DB
    % As an alternative, we could use simply coll.save(doc) without
    % writeconcern as the save method takes only one DBObject. This method is
    % generally used to update a document (retrieve it, modify it and save it).
    filtre_passe_bas_collection.insert(innerDoc1_Objects(j), wc1);
    pause(1);
    moyennage_collection.insert(innerDoc2_Objects(j), wc2);
    pause(1);
end

waitbar(100/100);
pause(0.1);
close(waitBarUploadFile);
    
set(hObject,'Enable', 'off');
% Update handles structure
guidata(hObject, handles);

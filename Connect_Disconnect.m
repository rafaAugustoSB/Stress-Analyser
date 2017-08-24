%% En tête

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DA SILVA BATISTA Rafael Augusto  %
% Download an upload of the data   %
% on the database (mongoDB)        %
% Last modification : 30/07/17     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Getting credencials 

[login, password] = logindlg(); % Getting login and password.

%% Connection to the database and downloading the data of the experience.

% in order to communicate to the DB (mongoDB) we need to use some commands
% from Java language, so we need to put the driver and import the library
% in our path.

% add the downloaded JAR file to java path
javaaddpath('mongo-java-driver-2.13.1.jar');

% Server info
myIp = '10.250.7.146';
myPort = 27017;
dataBaseName = 'dataECG';

% [login password] = logindlg();
% if isempty(login) || isempty(password)
%     return
% end

% import the library
import com.mongodb.*;

% Connect to remote server
m = Mongo(myIp, myPort);

% open DB
db = m.getDB(dataBaseName);

% Authentication
auth = db.authenticate(login, password);
if auth
    disp('Connection stablished ...');
else
    disp('Login or password incorrect ...');
    return
end

% get the list of collections
allCollecions = db.getCollectionNames(); % get collection name

% get a collection
ensemble_donnes_collection = db.getCollection('ensemble_donnees');

% finding the experience
dbObjTestEnsembleDonnees = BasicDBObject('nom_fichier_csv', 't1.csv');
myDoc_ensemble_donnes_collection = ensemble_donnes_collection.findOne(dbObjTestEnsembleDonnees);

dbObjTestCopyEnsembleDonnees = ensemble_donnes_collection.findOne(dbObjTestEnsembleDonnees).copy();

% get the the number of data per line
dataPerLigne = myDoc_ensemble_donnes_collection.get('frequence_utilisee');

% monitoring.
disp('En train de récupérer les données ...');

% get the indexof the objects for using in the other collection.
releves_donnees = myDoc_ensemble_donnes_collection.get('id_releves_donnees').toArray;
releves_donnees = cell2mat(cell(releves_donnees))';

% get a collection
releve_donnees_collection = db.getCollection('releve_donnees');

% Verifying if the data is already processed.
% get a collection
filtre_passe_bas_collection = db.getCollection('aaa_filtre_passe_bas');
% get a collection
moyennage_Coll = db.getCollection('aab_moyennage');

dbObjTestIfExist = BasicDBObject('id_', releves_donnees(1));

iddbObjTestIfExist = filtre_passe_bas_collection.findOne(dbObjTestIfExist);
if  isempty(iddbObjTestIfExist)
    disp('Cool! Le traitement va commencer ...');
else
    disp('Pas cool! Les données sont déjà traitées ...');
    return
end

% allocation of vectors
donnees_brutes = []; % vector of raw data.
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
        donnees_brutes = horzcat(donnees_brutes, donnesBrute_par_ligne);
    end
    
    % get the time of the start and end of the experience.
    if i == 1
        heure_debut = datetime(liste_dataTime(4, 1));
    end
    
    if i == length(releves_donnees)
        heure_fin = datetime(liste_dataTime(length(liste_dataTime), 1));
    end
    
    % monitoring
    disp('En train de récupérer les données ...');
end

% monitoring
disp('Données recuperées correctement ...');


%% PROCESSING ...

% Here we can call a file "traitement.m"

%% Upload de donnees vers la BDD.

% fist package of the signal processing - FILTRE PASSE BAS.

% vector of DBObjects.
innerDoc1_Objects = {};

% there is a limite concerning the size of the file upload, it is :
% 16793600 bit. So just like for uploading the file .csv in the DB we need to
% split the data too. It's is usefull because it keeps the same structure
% of other collections too (see releve_donnees).

data = 1;
k = 0;
for i = 1:length(qqte_donnees_par_objet)
    innerDoc1 = BasicDBObject();
    innerDoc1.put('id_', int32(releves_donnees(i)));
    innerDoc1.put('id_ensemble_donnees', int32(dbObjTestCopyEnsembleDonnees.get('id_')));
    for j = 1 : qqte_donnees_par_objet(i)
        innerDoc1.put(liste_date_complet(k + j,:), dataLowFiltre(data:(data + dataPerLigne - 1)));
        data = data + dataPerLigne -1;
    end
    k = k + qqte_donnees_par_objet(i);
    innerDoc1_Objects = [innerDoc1_Objects innerDoc1];
end

% second package of the signal processing - MOYENNAGE DATABASE.

% vector of DBObjects.
innerDoc2_Objects = {};

data = 1;
k = 0;
for i = 1:length(qqte_donnees_par_objet)
    innerDoc2 = BasicDBObject();
    innerDoc2.put('id_', int32(releves_donnees(i)));
    innerDoc2.put('id_ensemble_donnees', int32(dbObjTestCopyEnsembleDonnees.get('id_')));
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
    moyennage_collection.insert(innerDoc2_Objects(j), wc2);
end


%% Close the connection with the DB.

m.close;
m.close();
close(m);
disp('Connexion vers la BDD fermée ...');
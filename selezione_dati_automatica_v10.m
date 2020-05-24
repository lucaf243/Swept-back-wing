% Questo script apre in sequenza tutti i file esportati da XFLR5 e analizza
% i principali parametri aerodinamici al variare dell'angolo di freccia
% per un determinato angolo di attacco.
%
% Nell'analisi in XFRL5 sono stati tenuti costanti i seguenti parametri:
% - Velocità pari a 55 m/s
% - Angolo di attacco
% - Superficie alare, corda alare, AR, TR

result = isfolder('data_offset') ; 

a=any(size(dir(['data_offset' '/*.txt' ]),1));

if result==0
mkdir prova_data_offset
movefile MainWing_a=0.00_v=55.00ms_offset* data_offset
elseif result==1 && a==0
    movefile MainWing_a=0.00_v=55.00ms_offset* data_offset
end

clear a ans result

listafile = dir ('data_offset/*.txt');
[n_file,b] = size(listafile);

%calcolo degli offset
offset_3 = 0:0.5:5;
offset = offset_3;

%calcolo degli angoli di freccia


%inizializzazione dei vettori intensita del vortice 
 gamma_Tip_dx = zeros(n_file,1);
 gamma_Tip_sx = zeros(n_file,1);
 
 %inizializzazione del vettore velocità
 Vel = zeros(n_file,1);
 
 %inizializzo il vettore CP
 CP = zeros(n_file,3);
 
Cd=zeros; 
Cm=zeros; 
CL=zeros;

for i = 1:n_file
    
    %-------------CALCOLO INTENSITA DEI VORTICI ALLE ESTREMITA'------%
     filename = strcat("data_offset","/",listafile(i).name);
    delimiterIn    = ' ';
    headerlinesIn  = 21;
    format long;
    
    analysis_struct = importdata(filename,delimiterIn,headerlinesIn);
    
    
    %estraggo la velocita
    C = analysis_struct.textdata{4,1}; %estraggo la 4a riga
    C = C(1,10:15); %estraggo il valore della velocità
    Vel(i) = str2double(C);
    
    

    
    
    fields = fieldnames(analysis_struct);
    analysis  = char(fields(1));

    wing_analysis = analysis_struct.(analysis);
    
    [WA_row,WA_column] = size(wing_analysis);
    
    chord_Tip_dx = wing_analysis(1,2);
    Cl_WingTip_dx = wing_analysis(1,4);
    
    chord_Tip_sx = wing_analysis(WA_row,2);
    Cl_WingTip_sx = wing_analysis(WA_row,4);
    
    gamma_Tip_dx(i) = (Vel(i)*chord_Tip_dx*Cl_WingTip_dx)/2;
    gamma_Tip_sx(i) = (Vel(i)*chord_Tip_sx*Cl_WingTip_sx)/2;
    
    %----------------CATTURO GLI ALTRI DATI--------------------------------%
    
    %estraggo CL 
    C = analysis_struct.textdata{9,1}; 
    C = C(1,10:19); 
    CL(i) = str2double(C);
    
    %estraggo Cd 
    C = analysis_struct.textdata{11,1}; 
    C = C(1,10:19); 
    Cd(i) = str2double(C);
    
    %estraggo il Cm
    C = analysis_struct.textdata{13,1}; 
    C = C(1,8:18); 
    Cm(i) = str2double(C);
    
    %estraggo il CP
    C = analysis_struct.textdata{15,1}; 
    Cx = C(1,10:19); 
    Cy = C(1,34:43); 
    Cz = C(1,58:67); 
    CP(i,:) =[str2double(Cx) str2double(Cy) str2double(Cz)];

end


%prelevo semiapertura alare
Wing_Span = abs(wing_analysis(1,1));

%angoli di freccia relativi all'offset
tg_beta =atan(offset/Wing_Span);
tg_beta_deg = tg_beta*180/pi;

 
 E=CL./Cd; % efficienza

%calcolo fattore di osvald al variare dell'angolo di freccia
e = (CL.^2)./(pi*7.483.*Cd);                    % 7.483 è il valore dell'AR

% costruisco una matrice con tutti i coefficienti estratti da XFLR5
data_import=[tg_beta_deg' , CL' , Cd' , Cm' , gamma_Tip_dx , e']; 

%percentage ganches
ref = find(tg_beta_deg== 0);
CL_0 = CL(ref);
Cd_0 = Cd(ref);
gamma_Tip_dx_0 = gamma_Tip_dx(ref);
Cm_0 = Cm(ref);
e_0 = e(ref);

CL_variazioni = ((CL-CL_0)/CL_0)*100;
Cd_variazioni = ((Cd-Cd_0)/Cd_0)*100;
Cm_variazioni = ((Cm-Cm_0)/Cm_0)*100;
gamma_Tip_dx_variazioni = ((gamma_Tip_dx-gamma_Tip_dx_0)/gamma_Tip_dx_0)*100;
e_variazioni = ((e-e_0)/e_0)*100;

Tab_variazioni=[tg_beta_deg', CL_variazioni', Cd_variazioni', Cm_variazioni', ...
                    gamma_Tip_dx_variazioni, e_variazioni'];
% converto Tab_variazioni in formato TABLE (per esportarlo in un formato
% .txt,.xls,.csv)
T=array2table(Tab_variazioni,...
             'VariableNames',{'tg_beta_deg' 'CL_variazioni' 'Cd_variazioni' 'Cm_variazioni' 'gamma_tip_variaz' 'e_var'});

% costruisco una tabella con le variazioni minime e massime rispetto
% all'ala rettangolare

[r,c]=size(Tab_variazioni);                             % 1) estraggo la dimensione dell'array Tab_variazioni per il ciclo FOR

[var_min,var_max]=bounds(Tab_variazioni);               % 2) estraggo il valore massimo e minimo di ogni colonna di Tab_variazioni

Variazioni_max_min=zeros;                               % 3) inizializzo un array di zeri


for i=2:c
    
    [row1,~]=find(Tab_variazioni(:,i)==var_min(1,i));                % 4) cerco la riga corrispondente al valore minimo della colonna i 
                                                                        %
    Variazioni_max_min(2*(i-1),1)=var_min(i);                           % 5) associo alle righe dell'array Variazioni_max_min pari il valore minimo di ogni colonna
                                                                        %
    Variazioni_max_min(2*(i-1),2)=Tab_variazioni(row1,1);               % 6) estraggo il valore dell'angolo di freccia associato al valore minimo

    [row2,~]=find(Tab_variazioni(:,i)==var_max(1,i));                % procedimento analogo alla ricerca
     Variazioni_max_min(2*(i-1)+1,1)=var_max(i);                        % della variazione minima e dell'angolo 
     Variazioni_max_min(2*(i-1)+1,2)=Tab_variazioni(row2,1);            % di freccia corrispondente
     
end
Variazioni_max_min=Variazioni_max_min(2:end,:);                         % elimino la prima riga di zeri (bug del ciclo for)

% creo una tabella esportabile in altri formati

Variazioni_massime_minime=array2table(Variazioni_max_min, ...
                            'VariableNames', {'valore' 'angolo_di_freccia'}, ...
                            'RowNames', {'CL_min%' 'CL_max%' 'Cd_min%' 'Cd_max%' 'Cm_min%' 'Cm_max%' 'Gamma_min%' 'Gamma_max%' 'e_min%' 'e_max%'});

% per esportare la tabella aggiungere una delle seguenti righe:
% per formati .txt, .dat, .csv
% writetable(Variazioni_massime_minime,'Variazioni_massime_minime.dat','Delimiter',' ')
% per formati .xls, .xlsm, .xlsx
% writetable(Variazioni_massime_minime,'Variazioni_massime_minime.xls')

% valori desiderati CL_max ; Cd_min ; Cm_min ; Gamma_min ; e_max
% estraggo i valori di angolo di freccia associati ai suddetti valori

angoli_di_freccia=[Variazioni_max_min(2,2) ; Variazioni_max_min(3,2) ; Variazioni_max_min(5,2) ; ...
                    Variazioni_max_min(7,2) ; Variazioni_max_min(10,2)];

angoli_di_freccia=unique(angoli_di_freccia);                                     % controllo e rimuovo elementi uguali nel vettore angoli_di_freccia              
                
coefficienti=zeros(length(angoli_di_freccia),length(data_import(1,:)));        % inizializzo un vettore vuoto di coefficienti


for i=1:length(angoli_di_freccia)                                           %
                                                                            % coefficienti ha la seguente struttura
    [row,~]=find(data_import(:,1)==angoli_di_freccia(i));                   % coefficienti=[angoli_di_freccia, CL , Cd , E , Gamma , e]
                                                                            %
    coefficienti(i,1)=angoli_di_freccia(i);                                 %
    coefficienti(i,2:end)=data_import(row,2:end);                           %
    
end

%------------ASSEGNAZIONE-PUNTEGGI--------------------------%

punteggio(:,1)=coefficienti(:,1);   % prima colonna valore dell'angolo di freccia
                                    

    punteggio(:,2)=coefficienti(:,2)./(max(coefficienti(:,2)));                      % seconda colonna, assegno un punteggio per l'angolo di freccia che massimizza il CL
    punteggio(:,4)=abs(1-(coefficienti(:,4)./(max(coefficienti(:,4)))));             % quarta colonna, assegno un punteggio per l'angolo di freccia che minimizza il coefficiente di momento Cm
    punteggio(:,6)=coefficienti(:,6)./(max(coefficienti(:,6)));                      % sesta colonna, assegno un punteggio per l'angolo di freccia che massimizza il fattore di Oswald
    
    punteggio(:,3)=1-(coefficienti(:,3)./(max(coefficienti(:,3))));                  % terza colonna, assegno un punteggio per l'angolo di freccia che minimizza il Cd
    punteggio(:,5)=1-(coefficienti(:,5)./(max(coefficienti(:,5))));                  % quinta colonna, assegno un punteggio per l'angolo di freccia che minimizza l'intensita Gamma del vortice di estremità


punteggio1=punteggio(:,2:end);      % elimino la colonna contenente il valore dell'angolo di freccia
[row,~]=size(punteggio1);           % per sommare il punteggio totale ottenuto in ogni configurazione

somma_punteggi=zeros;               % inizializzo un vettore vuoto somma_punteggi

for i=1:row
    somma_punteggi(i)=sum(punteggio1(i,:));         % sommo i punteggi ottenuti
end
% caso in cui il vettore somma_punteggi ammette massimo con molteplicità
% pari ad 1
[~,row]=max(somma_punteggi);
disp(['l''angolo di freccia ottimale è pari a   ', num2str(punteggio(row,1)) ,'º']) 

% stampo i valori ottenuti per la configurazione ottimale
valori_config_ottimale=array2table(coefficienti(row,2:end),...
                        'VariableNames',{'CL' 'Cd_indotto' 'Cm' 'Gamma' 'e'});
disp(valori_config_ottimale)

% caso in cui il vettore somma_punteggi ammette un massimo con molteplicità
% maggiore di 1
if length(somma_punteggi)~=length(unique(somma_punteggi))
    warning('ci sono piu'' valori uguali nel vettore somma_punteggi, verificare se il valore con molteplicita>1 corrisponde al massimo')
end

punteggio_data_import(:,2)=data_import(:,2)./(max(data_import(:,2)));
punteggio_data_import(:,3)=1-(data_import(:,3)./(max(data_import(:,3))));
punteggio_data_import(:,4)=abs(1-(data_import(:,4)./(max(data_import(:,4)))));
punteggio_data_import(:,5)=1-(data_import(:,5)./(max(data_import(:,5))));
punteggio_data_import(:,6)=data_import(:,6)./(max(data_import(:,6)));

somma_punteggi_data_import=sum(punteggio_data_import,2);
[~,IDmax]=max(somma_punteggi_data_import);
disp(['l''angolo di freccia ottimale considerando tutti gli offset è pari a   ', num2str(data_import(IDmax,1)) ,'º'])

if data_import(IDmax,1)~=punteggio(row,1)

    warning('le due configurazioni ottimali non corrispondono')
end
%-----------FINE-ASSEGNAZIONE-PUNTEGGI--------------------------%

coefficienti1=coefficienti(:,2:end);

rownames = compose('%3.2f', angoli_di_freccia);

coefficienti_tab=array2table(coefficienti1, ...                                          % Coefficienti_tab è una tabella con il valori selezionati di
                            'VariableNames', {'CL' 'Cd' 'Cm' 'Gamma_tip' 'e'}, ...       % angolo di freccia ed i rispettivi CL,Cd,E,Gamma,e
                            'RowNames', rownames);                                       %

% elimino variabili inutili
 clear col col1 col2 i row row1 row2 row3 row4 row5   

%------------GRAFICI-----------------------------%
figure(1);clf;
plot(tg_beta_deg,gamma_Tip_dx,'r','LineWidth',2), hold on, plot(tg_beta_deg,gamma_Tip_sx,'b*')
grid on
title('intensità dei vortici di estremità alari al variare dell''angiolo di freccia')
legend('vortice destro','vortice sinistro')
xlabel('angolo di freccia (deg)')
ylabel('intensità del vortice')
hold off

figure(2);clf;
subplot(2,1,1)
hold on
grid on
plot(tg_beta_deg,CL,'LineWidth',2)
hold off
title('CL(offset)');
xlabel('angolo di freccia (deg)')

subplot(2,1,2)
hold on
grid on
plot(tg_beta_deg,Cd,'r','LineWidth',2,'color', [0.6350 , 0.0780 , 0.1840])
hold off
title('Cd(offset)');
xlabel('angolo di freccia (deg)')


figure(3);clf;
hold on
grid on
plot(tg_beta_deg,CL./Cd,'LineWidth',2)
hold off
title('efficienza(offset)');
xlabel('angolo di freccia (deg)')
ylabel('efficienza')


figure(4);clf;
hold on
grid on
plot(tg_beta_deg,e,'LineWidth',2)
hold off
title('fattore di Oswald (offset)')
xlabel('angolo di freccia (deg)')
ylabel('fattore di Oswald')

figure(5);clf;
hold on
grid on
plot(tg_beta_deg,Cm,'LineWidth',2,'color',[0, 0.5, 0])
hold off
title('coefficiente di momento')
xlabel('angolo di freccia(deg)')
ylabel('coefficiente di momento')

figure(6)
hold on
grid on
plot(tg_beta_deg,Cd,'LineWidth',2)
hold off
title('coefficiente di resistenza')
xlabel('angolo di freccia(deg)')
ylabel('Cd')
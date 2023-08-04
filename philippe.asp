% valeur absolue maximale prise par les variables = granularité qualitative
% le nombre de points de temps est pris égal à max+1 : -1 représente la période basale, 0 à max-1 la période d'intervention
 #const max = 21.

% valeur du point de temps où on décide de la continuation
 #const mean = max-1.

% valeur absolue minimale de la discontinuité de EI, dES et P entre -1 et 0
 #const disc = 1.

% bornes de l'écart positif entre E0 et Emin
 gapm(max/5).
 gapp(max/3).


% les cinq variables causales possibles
 cause((ei;des;es;ee;p)).

% signe des causalités
 1{csgn(I,(-1;0;1))}1 :- cause(I).

% les trois scénarios et leur existence temporelle (on garde la plage jusqu'à max-1 même si le scénario 
% finit avant car on veut que les conditions terminales sur les profils des variables s'appliquent toujours)
% ni ne dure que jusqu'à 0 : non initiation
% id dure jusqu'à mean : initiation puis discontinuation
% ind dure jusqu'à mean : initiation puis non discontinuation
 time(ni,-1..0).
 time(id,-1..max-1).
 time(ind,-1..max-1).


% profil de EI > 0 avec discontinuité entre -1 et 0 d'au moins disc unités et valeur finale moins que basale
% croissante et continue (variation de 0 ou 1 d'un point de temps au suivant) de 0 à max-2
% EI ne peut rester constante de 0 à mean
% 2*V1 < V0 :- v(Ph,-1,ei,V0); v(Ph,0,ei,V1).
% 2*(VM-V1) < (V0-V1) :- v(Ph,-1,ei,V0); v(Ph,0,ei,V1); v(Ph,max-1,ei,VM).
 1 { v(Ph,T,ei,1..max) } 1 :- time(Ph,T).
 V1 <= V0-disc :- v(Ph,-1,ei,V0); v(Ph,0,ei,V1).
 VM < V0 :- v(Ph,-1,ei,V0); v(Ph,max-1,ei,VM).
 VT <= VTP  :- v(Ph,T,ei,VT); v(Ph,TP,ei,VTP); time (Ph,T); T >= 0; T <= max-2; TP = T + 1.
 VTP <= VT+1 :- v(Ph,T,ei,VT); v(Ph,TP,ei,VTP); time (Ph,T); T >= 0; T <= max-2; TP = T + 1.
 Vm > V1 :- v(Ph,0,ei,V1); v(Ph,mean,ei,Vm).


% profil de dES ≤ 0, vaut 0 en -1 et à la fin, discontinuité d'au moins disc unités entre -1 et 0, 
% croissante et continue (variation de 0 ou 1 d'un point de temps au suivant) de 0 à max-2
% 2*(Vm-V1) > (V0-V1) :- v(Ph,-1,des,V0); v(Ph,0,des,V1); v(Ph,mean,des,Vm).
 v(Ph,-1,des,0) :- time(Ph,-1).
 1 { v(Ph,0,des,-max..-disc) } 1 :- time(Ph,0).
 1 { v(Ph,T,des,-max..0) } 1 :- time(Ph,T); T >= 1.
 v(Ph,max-1,des,0) :- time(Ph,max-1). 
 VT <= VTP  :- v(Ph,T,des,VT); v(Ph,TP,des,VTP); time (Ph,T); T >= 0; T <= max-2; TP = T + 1.
 VTP <= VT+1 :- v(Ph,T,des,VT); v(Ph,TP,des,VTP); time (Ph,T); T >= 0; T <= max-2; TP = T + 1.


% profil de ES > 0 continue décroissante (variation de 0 ou -1 d'un point de temps au suivant de -1 à max-2)
% ES ne peut rester constante de 0 à mean
% 2*(V0-VM) < V0 :- v(Ph,-1,es,V0); v(Ph,max-1,es,VM).
 1 { v(Ph,T,es,1..max) } 1 :- time(Ph,T).
 VT-1 <= VTP :- v(Ph,T,es,VT); v(Ph,TP,es,VTP); time (Ph,T); T <= max-2; TP = T + 1.
 VTP <= VT :- v(Ph,T,es,VT); v(Ph,TP,es,VTP); time (Ph,T); T <= max-2; TP = T + 1.
 Vm < V1 :- v(Ph,0,es,V1); v(Ph,mean,es,Vm).
% si dES vaut 0 à un point de temps, ES ne change pas de valeur entre ce point de temps et le suivant
  V = V1 :- v(Ph,T,des,0); v(Ph,T,es,V); v(Ph,T+1,es,V1).


% profil de EE > 0 continue décroissante (variation de 0 ou -1 d'un point de temps au suivant de -1 à max-2, pas de variation de -1 à 0)
% EE ne peut rester constante de 0 à mean
% 2*(V0-VM) < V0 :- v(Ph,-1,ee,V0); v(Ph,max-1,ee,VM).
 1 { v(Ph,T,ee,1..max) } 1 :- time(Ph,T).
 V0 = V1 :- v(Ph,-1,ee,V0); v(Ph,0,ee,V1).
 VT-1 <= VTP  :- v(Ph,T,ee,VT); v(Ph,TP,ee,VTP); time (Ph,T); T >= 0; T <= max-2; TP = T + 1.
 VTP <= VT :- v(Ph,T,ee,VT); v(Ph,TP,ee,VTP); time (Ph,T); T >= 0; T <= max-2; TP = T + 1.
 Vm < V1 :- v(Ph,0,ee,V1); v(Ph,mean,ee,Vm).


% profil de P, 0 en -1 puis constante > 0 de 0 à max-1 avec discontinuité d'au moins disc unités entre -1 et 0
 v(Ph,-1,p,0) :- time(Ph,-1).
 v(Ph,T,p,V) :- time(Ph,T), T >= 0, val(p,Ph,V).
 1 { val(p,Ph, disc..max) } 1 :- time(Ph,_).


% Application du signe de la causalité
 eff(Ph,T,I,V*W) :- v(Ph,T,I,V); csgn(I,W).

% les trois scénarios
% ni : non initiation à 0
 :- gapp(G); not #sum { W,I,0: eff(ni,0,I,W); -W,I,-1: eff(ni,-1,I,W)  } -G.

% id : initiation à 0 suivie d'une discontinuation entre 1 et mean
  :-  gapm(G); not 1-G #sum { W,I,0: eff(id,0,I,W); -W,I,-1: eff(id,-1,I,W)  }.
  :- not 1 #count {T : 1 <= T, T <= mean - 1, dis(T) }.
   dis(T) :- gapp(G); time(id,T), #sum { W,I,T: eff(id,T,I,W); -W,I,-1: eff(id,-1,I,W)  } -G.

% ind : initiation à 0 suivie de non discontinuation jusqu'à mean
  :- time(ind,T); gapm(G); T <= mean - 1, not 1-G #sum { W,I,T+1: eff(ind,T+1,I,W); -W,I,-1: eff(ind,-1,I,W) }.
% effet cumulé sur E est le même à la fin qu'en basal 
  :- not  #sum { W,I,max-1: eff(ind,max-1,I,W); -W,I,-1: eff(ind,-1,I,W)} = 0.


% Contrainte : l'effet cumulé sur E des causalités est > gapp en -1
 :- time(Ph,_); gapp(G); not G+1 #sum { W,I: eff(Ph,-1,I,W) }.


% au plus 3 causalités actives sur les 5 possibles
 :- not 2 #count { I : csgn(I,0) }.


% On exclut d'avoir trois causalités directes actives de EI, dES, EE 
% vu que les 3 variables sont liées par une équation (égalité) de bilan thermodynamique.
% :- csgn(ei,Wei); csgn(des,Wdes); csgn(ee,Wee); Wei != 0; Wdes != 0; Wee != 0.

% On exclut d'avoir deux causalités directes actives de dES et ES 
% vu que les deux variables sont liées par la dérivée
% :- csgn(des,Wdes); csgn(es,Wes); Wdes != 0; Wes != 0.

% On exclut d'avoir deux causalités directes actives de ES, EE.
% :- csgn(es,Wes); csgn(ee,Wee); Wes != 0; Wee != 0.

% On exclut d'avoir deux causalités directes actives de dES et EE 
% :-csgn(des,Wdes); csgn(ee,Wee); Wdes != 0; Wee != 0.


 #show csgn/2.

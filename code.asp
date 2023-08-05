% maximal absolute value of the variables
	#const m = 21.
% value of the timepoint where non-discontinuation or continuation are decided 
 	#const mean = m-1.
% minimal absolute value of the discontinuity of EI, dES, and P
 	#const disc = 1.
% bounds of the positive gap between E0 and Emin
 	deltaem(m/5).
 	deltaep(m/3).

% the five possible causal input variables to psychological energy
 	cause((ei;des;es;ee;p)).
% sign of the causalities
 	1{csgn(I,(-1;0;1))}1 :- cause(I).
% the three core phenomena and their temporal extent
% number of timepoints is m+1: -1 for the basal period, 0 to m-1 for the intervention period
% ni ends at time 0, id and ind at time mean, but m-1 is used because the conditions on the 
% evolution profiles of the variables have to be taken into account up to the end
 	time(ni,-1..0).
 	time(id,-1..m-1).
 	time(ind,-1..m-1).

% evolution profile of EI > 0 with discontinuity of at least disc units between -1 and 0,
% then continuously increasing from 0 to the end (variation by 0 or 1 from a timepoint to the following),
% of at least one unit when reaching mean, without reaching its basal value at the end
 	1 { v(Ph,T,ei,1..m) } 1 :- time(Ph,T).
 	V1 <= V0-disc :- v(Ph,-1,ei,V0); v(Ph,0,ei,V1).
 	VM < V0 :- v(Ph,-1,ei,V0); v(Ph,m-1,ei,VM).
 	VT <= VTP  :- v(Ph,T,ei,VT); v(Ph,TP,ei,VTP); time (Ph,T); T >= 0; T <= m-2; TP = T + 1.
 	VTP <= VT+1 :- v(Ph,T,ei,VT); v(Ph,TP,ei,VTP); time (Ph,T); T >= 0; T <= m-2; TP = T + 1.
 	Vm > V1 :- v(Ph,0,ei,V1); v(Ph,mean,ei,Vm).
% evolution profile of dES <= 0, null in -1 and at the end, with discontinuity of at least disc units between -1 and 0,
% continuously increasing from 0 to the end (variation by 0 or 1 from a timepoint to the following)
 	v(Ph,-1,des,0) :- time(Ph,-1).
 	1 { v(Ph,0,des,-m..-disc) } 1 :- time(Ph,0).
 	1 { v(Ph,T,des,-m..0) } 1 :- time(Ph,T); T >= 1.
 	v(Ph,m-1,des,0) :- time(Ph,m-1). 
 	VT <= VTP  :- v(Ph,T,des,VT); v(Ph,TP,des,VTP); time (Ph,T); T >= 0; T <= m-2; TP = T + 1.
 	VTP <= VT+1 :- v(Ph,T,des,VT); v(Ph,TP,des,VTP); time (Ph,T); T >= 0; T <= m-2; TP = T + 1.
% evolution profile of ES > 0 continuously decreasing from -1 to the end (variation by 0 or -1 
% from a timepoint to the following), of at least one unit between 0 and mean
 	1 { v(Ph,T,es,1..m) } 1 :- time(Ph,T).
 	VT-1 <= VTP :- v(Ph,T,es,VT); v(Ph,TP,es,VTP); time (Ph,T); T <= m-2; TP = T + 1.
 	VTP <= VT :- v(Ph,T,es,VT); v(Ph,TP,es,VTP); time (Ph,T); T <= m-2; TP = T + 1.
 	Vm < V1 :- v(Ph,0,es,V1); v(Ph,mean,es,Vm).
% if dES is null at a timepoint, ES will not change value between this timepoint and the following
  	V = V1 :- v(Ph,T,des,0); v(Ph,T,es,V); v(Ph,T+1,es,V1).
% evolution profile of EE > 0 continuously decreasing from -1 to the end (variation by 0 or -1 
% from a timepoint to the following), of at least one unit between 0 and mean
 	1 { v(Ph,T,ee,1..m) } 1 :- time(Ph,T).
 	V0 = V1 :- v(Ph,-1,ee,V0); v(Ph,0,ee,V1).
 	VT-1 <= VTP  :- v(Ph,T,ee,VT); v(Ph,TP,ee,VTP); time (Ph,T); T >= 0; T <= m-2; TP = T + 1.
 	VTP <= VT :- v(Ph,T,ee,VT); v(Ph,TP,ee,VTP); time (Ph,T); T >= 0; T <= m-2; TP = T + 1.
 	Vm < V1 :- v(Ph,0,ee,V1); v(Ph,mean,ee,Vm).
% profil de P, 0 en -1 puis constante > 0 de 0 à m-1 avec discontinuité d'au moins disc unités entre -1 et 0
% evolution profile of P >= 0, null in -1, with discontinuity 
% of at least disc units between -1 and 0, then constant 
 	v(Ph,-1,p,0) :- time(Ph,-1).
 	v(Ph,T,p,V) :- time(Ph,T), T >= 0, val(p,Ph,V).
 	1 { val(p,Ph, disc..m) } 1 :- time(Ph,_).

% application of the sign of the causality to get the effect
 	eff(Ph,T,I,V*W) :- v(Ph,T,I,V); csgn(I,W).

% the three core phenomena
% Ph.1: non-initiation (ni)
 	:- deltaep(G); not #sum { W,I,0: eff(ni,0,I,W); -W,I,-1: eff(ni,-1,I,W)  } -G.
% Ph.2: initiation followed by discontinuation (id) between 1 and mean
  	:- deltaem(G); not 1-G #sum { W,I,0: eff(id,0,I,W); -W,I,-1: eff(id,-1,I,W)  }.
  	:- not 1 #count {T : 1 <= T, T <= mean - 1, dis(T) }.
  	dis(T) :- deltaep(G); time(id,T), #sum { W,I,T: eff(id,T,I,W); -W,I,-1: eff(id,-1,I,W)  } -G.
% Ph.3: initiation followed by non-discontinuation (ind) up to mean
  	:- time(ind,T); deltaem(G); T <= mean - 1, not 1-G #sum { W,I,T+1: eff(ind,T+1,I,W); -W,I,-1: eff(ind,-1,I,W) }.
% the global causal effect on E at -1 is greater than deltaep
 	:- time(Ph,_); deltaep(G); not G+1 #sum { W,I: eff(Ph,-1,I,W) }.

% restoration at the end of the basal value of E
  	:- not #sum { W,I,m-1: eff(ind,m-1,I,W); -W,I,-1: eff(ind,-1,I,W)} = 0.

% minimality: at most three active causalities (change 2 by 3 if restoration disabled)
 :- not 2 #count { I : csgn(I,0) }.

% display of the signs of the causalities
 	#show csgn/2.

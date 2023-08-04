% maximal absolute value of the variables
	#const m=21.
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
% (-1 for basal, 0 for initiation, 1 for discontinuation periods)
 	time(ni,-1..0).
 	time(id,-1..1).
 	time(ind,-1..1).

% evolution profile of EI > 0 with discontinuity of at least disc  
% units between -1 and 0 then rise of at least one unit without 
% reaching its basal value
 	1 { v(Ph,T,ei,1..m) } 1 :- time(Ph,T).
 	V2 <= V1-disc :- v(Ph,-1,ei,V1); v(Ph,0,ei,V2).
 	V3 > V2 :- v(Ph,0,ei,V2); v(Ph,1,ei,V3).
 	V3 < V1 :- v(Ph,-1,ei,V1); v(Ph,1,ei,V3).
% evolution profile of dES <= 0 with discontinuity of at least disc  
% units between -1 and 0, null in -1 and in 1 
 	v(Ph,-1,des,0) :- time(Ph,-1).
 	1 { v(Ph,0,des,-m..-disc) } 1 :- time(Ph,0).
 	v(Ph,1,des,0) :- time(Ph,1). 
% evolution profile of ES > 0 decreasing by at least one unit with  
% a variation of at most one unit (continuity) between -1 and 0 
 	1 { v(Ph,T,es,1..m) } 1 :- time(Ph,T).
 	V1 >= V2  :- v(Ph,-1,es,V1); v(Ph,0,es,V2).
 	V2 >= V1-1  :- v(Ph,-1,es,V1); v(Ph,0,es,V2).
 	V1 > V3 :- v(Ph,-1,es,V1); v(Ph,1,es,V3).
% evolution profile of EE > 0 decreasing by at least one unit with  
% a variation of at most one unit (continuity) between -1 and 0 
 	1 { v(Ph,T,ee,1..m) } 1 :- time(Ph,T).
 	V1 >= V2  :- v(Ph,-1,ee,V1); v(Ph,0,ee,V2).
 	V2 >= V1-1   :- v(Ph,-1,ee,V1); v(Ph,0,ee,V2).
 	V1 > V3 :- v(Ph,-1,ee,V1); v(Ph,1,ee,V3).
% evolution profile of P >= 0, null in -1, with discontinuity 
% of at least disc units between -1 and 0 
 	v(Ph,-1,p,0) :- time(Ph,-1).
 	v(Ph,T,p,V) :- time(Ph,T); T >= 0; val(p,Ph,V).
 	1 { val(p,Ph, disc..m) } 1 :- time(Ph,_).

% application of the sign of the causality to get the effect
 	eff(Ph,T,I,V*W) :- v(Ph,T,I,V); csgn(I,W).

% the three core phenomena
% Ph.1: non-initiation (ni)
  	:- deltaep(G); 
   not #sum { W,I,0: eff(ni,0,I,W); -W,I,-1: eff(ni,-1,I,W) } -G.
% Ph.2: initiation followed by discontinuation (id)
  	:- deltaem(G); 
   not 1-G #sum { W,I,0: eff(id,0,I,W); -W,I,-1: eff(id,-1,I,W)}.
  	:- deltaep(G); 
   not #sum { W,I,1 : eff(id,1,I,W); -W,I,-1: eff(id,-1,I,W)} -G.
% Ph.3: initiation followed by non-discontinuation (ind)
  	:- deltaem(G); 
   not 1-G #sum { W,I,0: eff(ind,0,I,W); -W,I,-1: eff(ind,-1,I,W)}.
 	:- deltaem(G); 
   not 1-G #sum { W,I,1: eff(ind,1,I,W); -W,I,-1: eff(ind,-1,I,W)}.
% the global causal effect on E at -1 is greater than deltaep
 	:- time(Ph,_); deltaep(G); 
   not G+1 #sum { W,I: eff(Ph,-1,I,W) }.

% restoration at 1 of the basal value of E
  	:- not  
   #sum { W,I,1: eff(ind,1,I,W); -W,I,1: eff(ind,-1,I,W)} = 0.

% minimality: at most three active causalities
 	:- not 2 #count { I : csgn(I,0) }.
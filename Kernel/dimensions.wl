
(*   This file is part of UnitSystems                   *)
(*   It is licensed under the MIT license               *)
(*   UnitSystems Copyright (C) 2021 Michael Reed        *)
(*       _           _                         _        *)
(*      | |         | |                       | |       *)
(*   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _  *)
(*  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` | *)
(* | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| | *)
(*  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_| *)
(*                                                      *)
(*   https://github.com/chakravala                      *)
(*   https://crucialflow.com                            *)

usq = {"F", "M", "L","T","Q","\[Theta]","N","J","A","\[CapitalLambda]","C"}

d1 = USQ[{0,0,0,0,0,0,0,0,0,0,0},1]
dF = USQ[{1,0,0,0,0,0,0,0,0,0,0},1]
dM = USQ[{0,1,0,0,0,0,0,0,0,0,0},1]
dL = USQ[{0,0,1,0,0,0,0,0,0,0,0},1]
dT = USQ[{0,0,0,1,0,0,0,0,0,0,0},1]
dQ = USQ[{0,0,0,0,1,0,0,0,0,0,0},1]
d0 = USQ[{0,0,0,0,0,1,0,0,0,0,0},1]
dN = USQ[{0,0,0,0,0,0,1,0,0,0,0},1]
dJ = USQ[{0,0,0,0,0,0,0,1,0,0,0},1]
dA = USQ[{0,0,0,0,0,0,0,0,1,0,0},1]
dR = USQ[{0,0,0,0,0,0,0,0,0,1,0},1]
dC = USQ[{0,0,0,0,0,0,0,0,0,0,1},1]

usqbox[{x_, 1}] := ToBoxes[x]
usqbox[{x_, 1.}] := ToBoxes[x]
usqbox[{x_, y_}] := SuperscriptBox[x, ToBoxes[y]]
usqbox[{_, 0}] = Sequence[]
usqbox[{_, 0.}] = Sequence[]

USQ /: MakeBoxes[x_USQ, StandardForm] := unitBoxes[x]

USQ /: Plus[a:USQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_],USQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},y_]] := USQ[z,x+y]
USQ /: Subtract[USQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_],USQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},y_]] := USQ[z,x-y]
USQ /: Times[USQ[a_List,x_],USQ[b_List,y_]] := USQ[a+b,x*y]
USQ /: Divide[USQ[a_List,x_],USQ[b_List,y_]] := USQ[a-b,x/y]
USQ /: Power[USQ[a_List,c_],b_] := USQ[b*a,Power[c,b]]
USQ /: Sqrt[USQ[a_List,c_]] := USQ[a/2,Sqrt[c]]
USQ /: Inverse[USQ[a_List,c_]] := USQ[-a,1/c]

USQ /: Times[USQ[a_List,c_],b_] := USQ[a,c*b]

ToUSQ[Length] = dL
ToUSQ[Area] = dL^2
ToUSQ[Volume] = dL^3
ToUSQ[Power] = dF*dL/dT
ToUSQ[MomentOfInertia] = dM*dL^2
ToUSQ[SolidAngle] = dA^2

(*USQ /: Times[Length, a_USQ] := Times[ToUSQ[Length], a];
USQ /: Times[a_USQ, Length] := Times[a, ToUSQ[Length]];
USQ /: Divide[a_USQ, Length] := Divide[a, ToUSQ[Length]];
USQ /: Divide[Length, a_USQ] := Divide[ToUSQ[Length], a];

USQ /: Times[Area, a_USQ] := Times[ToUSQ[Area], a];
USQ /: Times[a_USQ, Area] := Times[a, ToUSQ[Area]];
USQ /: Divide[a_USQ, Area] := Divide[a, ToUSQ[Area]];
USQ /: Divide[Area, a_USQ] := Divide[ToUSQ[Area], a];

USQ /: Times[Volume, a_USQ] := Times[ToUSQ[Volume], a];
USQ /: Times[a_USQ, Volume] := Times[a, ToUSQ[Volume]];
USQ /: Divide[a_USQ, Volume] := Divide[a, ToUSQ[Volume]];
USQ /: Divide[Volume, a_USQ] := Divide[ToUSQ[Volume], a];

USQ /: Times[Power, a_USQ] := Times[ToUSQ[Power], a];
USQ /: Times[a_USQ, Power] := Times[a, ToUSQ[Power]];
USQ /: Divide[a_USQ, Power] := Divide[a, ToUSQ[Power]];
USQ /: Divide[Power, a_USQ] := Divide[ToUSQ[Power], a];

USQ /: Times[MomentOfInertia, a_USQ] := Times[ToUSQ[MomentOfInertia], a];
USQ /: Times[a_USQ, MomentOfInertia] := Times[a, ToUSQ[MomentOfInertia]];
USQ /: Divide[a_USQ, MomentOfInertia] := Divide[a, ToUSQ[MomentOfInertia]];
USQ /: Divide[MomentOfInertia, a_USQ] := Divide[ToUSQ[MomentOfInertia], a];

USQ /: Times[SolidAngle, a_USQ] := Times[ToUSQ[SolidAngle], a];
USQ /: Times[a_USQ, SolidAngle] := Times[a, ToUSQ[SolidAngle]];
USQ /: Divide[a_USQ, SolidAngle] := Divide[a, ToUSQ[SolidAngle]];
USQ /: Divide[SolidAngle, a_USQ] := Divide[ToUSQ[SolidAngle], a];*)

(* Measure *)

showGroup[u_] := If[AbstractUnitSystemQ[u],usq,showGroup[StringJoin["Abstract",u]]]
showGroup["AbstractUnified"] := {"kB","\[HBar]","c","\[Mu]0","me","Mu","Kcd","\[Phi]","\[Lambda]","\[Alpha]L","g0"}
showGroup["AbstractMetric"] := {"kgf", "kg", "m","s","C","K","mol","lm","rad","",""}
showGroup["AbstractMeridian"] := {"kegf", "keg", "em","s","eC","K","eg-mol","lm","rad","",""}
showGroup["AbstractBritish"] := {"lb", "slug", "ft","s","C","\[Degree]R","slug-mol","lm","rad","",""}
showGroup["AbstractEnglish"] := {"lbf", "lbm", "ft","s","C","\[Degree]R","lb-mol","lm","rad","",""}
showGroup["AbstractIPS"] := {"lb", "slinch", "in","s","C","\[Degree]R","slinch-mol","lm","rad","",""}
showGroup["AbstractFPS"] := {"pdl", "lb", "ft","s","C","\[Degree]R","lb-mol","lm","rad","",""}
showGroup["AbstractGauss"] := {"gf", "g", "cm","s","_","K","mol","lm","rad","",""}
showGroup["AbstractIAU"] := {"M\[CircleDot]f", "M\[CircleDot]", "au","D","C","K","mol","lm","rad","",""}
showGroup["AbstractIAUE"] := {"MEf", "ME", "au","D","C","K","mol","lm","rad","",""}
showGroup["AbstractIAUJ"] := {"MJf", "MJ", "au","D","C","K","mol","lm","rad","",""}
showGroup["AbstractMTS"] := {"tf", "t", "m","s","C","K","mol","lm","rad","",""}
showGroup["AbstractKKH"] := {"kgf", "kg", "km","h","C","K","mol","lm","rad","",""}
showGroup["AbstractMPH"] := {"lbf", "lb", "mi","h","C","\[Degree]R","lb-mol","lm","rad","",""}
showGroup["AbstractNautical"] := {"kegf", "keg", "nm","h","eC","K","mol","lm","rad","",""}
showGroup["AbstractFFF"] := {"firf", "fir", "fur","ftn","inf","\[Degree]R","fir-mol","lm","rad","",""}

showGroup["AbstractSI2019"] := showGroup["Metric"]
showGroup["AbstractSI1976"] := showGroup["Metric"]
showGroup["AbstractCODATA"] := showGroup["Metric"]
showGroup["AbstractConventional"] := showGroup["Metric"]
showGroup["AbstractInternational"] := showGroup["Metric"]
showGroup["AbstractInternationalMean"] := showGroup["Metric"]
showGroup["AbstractMetricEngineering"] := showGroup["Metric"]
showGroup["AbstractGravitationalMetric"] := showGroup["Metric"]
showGroup["AbstractSI2019Engineering"] := showGroup["MetricEngineering"]
showGroup["AbstractGravitationalSI2019"] := showGroup["GravitationalMetric"]
showGroup["AbstractMeridianEngineering"] := showGroup["Meridian"]
showGroup["AbstractGravitationalMeridian"] := showGroup["Meridian"]
showGroup["AbstractSurvey"] := showGroup["English"]
showGroup["AbstractEMU"] := showGroup["Gauss"]
showGroup["AbstractESU"] := showGroup["Gauss"]
showGroup["AbstractLorentzHeaviside"] := showGroup["Gauss"]
showGroup["AbstractKennelly"] := showGroup["Metric"]

unitBoxes[USQ[{0,0,0,0,0,0,0,0,0,0,0}, 1],_:usq] := "\[DoubleStruckL]"
unitBoxes[USQ[l_List, 1],s_:usq] := RowBox[Map[usqbox, Transpose[{s, l}]]]
unitBoxes[USQ[l_List, c_],s_:usq] :=  RowBox[Prepend[Map[usqbox, Transpose[{s, l}]], ToBoxes[c]]]

transformBoxes[d_,u_] := unitBoxes[transform[u,d],showGroup[u]]

Measure /: MakeBoxes[Measure[v_, d_USQ, u_String],StandardForm] := RowBox[{ToBoxes[v], "[", transformBoxes[d,u], "]",  u}]

Measure /: Times[Measure[v_,d_,u_],a_] := Measure[v*a,d,u]
Measure /: Plus[Measure[v_,d1,u_],a_] := Measure[v+a,d1,u]
Measure /: Subtract[Measure[v_,d1,u_],a_] := Measure[v-a,d1,u]
Measure /: Subtract[a_,Measure[v_,d1,u_]] := Measure[a-v,d1,u]

Measure /: Times[Measure[v1_,d1_,u_],Measure[v2_,d2_,u_]] := Measure[v1*v2,d1*d2,u]
Measure /: Divide[Measure[v1_,d1_,u_],Measure[v2_,d2_,u_]] := Measure[v1/v2,d1/d2,u]
Measure /: Power[Measure[v_,d_,u_],b_] := Measure[v^b,d^b,u]
Measure /: Plus[Measure[v1_,d_,u_],Measure[v2_,d_,u_]] := Measure[v1+v2,d,u]
Measure /: Subtract[Measure[v1_,d_,u_],Measure[v2_,d_,u_]] := Measure[v1-v2,d,u]
Measure /: Sqrt[Measure[v_,d_,u_]] := Measure[Sqrt[v],Sqrt[d],u]

Measure /: Equal[Measure[v_,USQ[l_,c_],u_],a_] := And[Total[l]==0,a==v*c]

(*quantity /: Log[b_,Power[B_,_USQ]] := "test"*)

expo[b_,e_] := If[e==0,1,b^e]
ratio[d_, u_, s_] := ratio[d,MatchSystem[s,u],UnitSystem[s]]
ratio[d_, u_UnitSystem, s_UnitSystem] := ratiocalc[TRANSFORM[d], u, s]
ratiocalc[USQ[d_,c_], u_, s_] := Unit[
	expo[BoltzmannConstant[u, s],d[[1]]]*
	expo[ReducedPlanckConstant[u, s],d[[2]]]*
	expo[SpeedOfLight[u, s],d[[3]]]*
	expo[MagneticConstant[u, s],d[[4]]]*
	expo[ElectronMass[u, s],d[[5]]]*
	expo[MonochromaticRadiation540THzLuminousEfficacy[u, s],d[[6]]]*
	expo[MolarMassConstant[u, s],d[[7]]]*
	expo[AngleConstant[u, s],d[[8]]]*
	expo[RationalizationConstant[u, s],d[[9]]]*
	expo[LorentzConstant[u, s],d[[10]]]*
	expo[GravityConstant[u, s],d[[11]]]]

(d_USQ)[u_String, s_String] := ratio[d,u,s]
(d_USQ)[u_UnitSystem,s_UnitSystem] := ratio[d,u,s]
(d_USQ)[v_,u_String] := d[v,UnitSystem[u]]
(d_USQ)[v_, u_String, s_String] := d[v, UnitSystem[u], UnitSystem[s]]
(d_USQ)[u_String] := d[UnitSystem[u]]
(*(d_USQ)[u_String, s_String] := d[UnitSystem[u], UnitSystem[s]]*)
(d_USQ)[v_, u_UnitSystem] := d[v, u, DefaultSystem[u]]
(d_USQ)[v_, u_UnitSystem, s_UnitSystem] := Module[{n = d[u, s]}, If[OneQ[n], v, v/n]]
(d_USQ)[u_UnitSystem] := d[u, DefaultSystem[u]]

Measure /: Entity["UnitSystem", s_][Measure[v_, d_, u_]] := Measure[v*ratio[d, u, s], d, s]
Measure[v_, d_, u_][Entity["UnitSystem",s_]] := Measure[v*ratio[d, u, s], d, s]
Measure[v_, d_, u_][s_String] := Measure[v*ratio[d, u, s], d, s]
ConvertUnit[Measure[v_, d_, u_],s_] := Measure[v*ratio[d, u, s], d, s]

dimension[Measure[_,d_,_]] := d

USQ[d_Symbol] := dimension[d[MM /@ UnitSystem["Natural"], UnitSystem["Natural"]]]

(* 1,2,3,4, 5, 6, 7,  8,9,10,11 *)
(*kB,ħ,𝘤,μ₀,mₑ,Mᵤ,Kcd,A,λ,αL,g₀ *)
(* F,M,L,T, Q, Θ, N,  J,A,Λ, C  *)

TRANSFORM[USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] := USQ[{-h,l+t+q/2-f-j,3*f+2*h+4*j-l-2*t-q/2,q/-2,m+h+n+2*(f+j)-l-t,-n,j,l+t+a+q/2-f-j,r-q/2,-q-c,l+t-h-2*(f+j)},x]

transform[u_String,x_] := If[AbstractUnitSystemQ[u],transform[StringTrim[u,"Abstract"],x],transform["Metric",x]]

transform["Unified",u_USQ] := TRANSFORM[u]
transform["MetricEngineering",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{f,m,l,t,q,h,n,j,a,0,0},x]
transform["GravitationalMetric",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{f+m,0,l-m,t+2*m,q,h,n,j,0,0,0},x]
transform["Metric",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,f+m,f+l,t-2*f,q,h,n,j,0,0,0},x]
(*transform["Astronomical",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{f,0,l+3*m,t-2*m,q,h,n,j,a,0,0},x]*)

transform["Gauss",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,f+m+q/2,f+l+(3/2)*q+c,t-2*(f+q)-c,0,h,n,j,0,0,0},x]
transform["ESU",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,f+m+q/2,f+l+(3/2)*q,t-2*(f+q),0,h,n,j,0,0,0},x]
transform["EMU",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,f+m+q/2,f+l+q/2,t-q-2*f,0,h,n,j,0,0,0},x]

transform["Stoney",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,f+m+h,0,l+t-f,q,0,0,0,0,0,0},x]
transform["Electronic",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,0,0,l+t-f,q,0,0,0,0,0,0},x]
transform["QCDoriginal",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,2*f+m+h-l-t,0,0,q,0,0,0,a,0,0},x]
transform["PlanckGauss",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,2*f+m+h-l-t,0,0,0,0,0,0,a,0,0},x]
transform["Planck",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,2*f+m+q+h-l-t,0,0,0,0,0,0,0,0,0},x]
transform["Natural",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,0,0,0,0,0,0,0,0,0,0},x]
transform["NaturalGauss",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,0,0,0,0,0,0,0,a,0,0},x]

transform["Rydberg",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,m+t+h-f,0,l+2*(t-h)-3*f,q,0,0,0,0,0,0},x]
transform["Hartree",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,0,l+2*(t-h)-3*f,0,q,0,0,0,0,0,0},x]
transform["Hubble",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,f+m,l+t-f,0,q,h,n,j,0,0,0},x]
transform["CosmologicalQuantum",USQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	USQ[{0,m+2*f-l-t,0,0,q,h,n,j,0,0,0},x]

transform["LorentzHeaviside",x_] := transform["Gauss",x]
(*transform["Thomson",x_] := transform["EMU",x]*)
transform["Kennelly",x_] := transform["EMU",x]
transform["Schrodinger",x_] := transform["Rydberg",x]
transform["QCD",x_] := transform["Planck",x]
transform["QCDGauss",x_] := transform["PlanckGauss",x]
transform["Cosmological",x_] := transform["Hubble",x]

transform["SI2019Engineering",x_] := transform["MetricEngineering",x]
transform["GravitationalSI2019",x_] := transform["GravitationalMetric",x]
transform["MeridianEngineering",x_] := transform["MetricEngineering",x]
transform["GravitationalMeridian",x_] := transform["GravitationalMetric",x]
transform["British",x_] := transform["GravitationalMetric",x]
transform["IPS",x_] := transform["GravitationalMetric",x]
transform["English",x_] := transform["MetricEngineering",x]
transform["Survey",x_] := transform["MetricEngineering",x]


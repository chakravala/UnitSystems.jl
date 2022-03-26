
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

isq = {"F", "M", "L","T","Q","\[Theta]","N","J","A","\[CapitalLambda]","C"}

d1 = ISQ[{0,0,0,0,0,0,0,0,0,0,0},1]
dF = ISQ[{1,0,0,0,0,0,0,0,0,0,0},1]
dM = ISQ[{0,1,0,0,0,0,0,0,0,0,0},1]
dL = ISQ[{0,0,1,0,0,0,0,0,0,0,0},1]
dT = ISQ[{0,0,0,1,0,0,0,0,0,0,0},1]
dQ = ISQ[{0,0,0,0,1,0,0,0,0,0,0},1]
d0 = ISQ[{0,0,0,0,0,1,0,0,0,0,0},1]
dN = ISQ[{0,0,0,0,0,0,1,0,0,0,0},1]
dJ = ISQ[{0,0,0,0,0,0,0,1,0,0,0},1]
dA = ISQ[{0,0,0,0,0,0,0,0,1,0,0},1]
dR = ISQ[{0,0,0,0,0,0,0,0,0,1,0},1]
dC = ISQ[{0,0,0,0,0,0,0,0,0,0,1},1]

isqbox[{x_, 1}] := ToBoxes[x]
isqbox[{x_, 1.}] := ToBoxes[x]
isqbox[{x_, y_}] := SuperscriptBox[x, ToBoxes[y]]
isqbox[{_, 0}] = Sequence[]
isqbox[{_, 0.}] = Sequence[]

ISQ /: MakeBoxes[ISQ[{0,0,0,0,0,0,0,0,0,0,0}, 1], StandardForm] := "\[DoubleStruckL]"
ISQ /: MakeBoxes[ISQ[l_List, 1], StandardForm] := RowBox[Map[isqbox, Transpose[{isq, l}]]]
ISQ /: MakeBoxes[ISQ[l_List, c_], StandardForm] :=  RowBox[Prepend[Map[isqbox, Transpose[{isq, l}]], c]]

ISQ /: Plus[a:ISQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_],ISQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},y_]] := ISQ[z,x+y]
ISQ /: Subtract[ISQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_],ISQ[z:{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},y_]] := ISQ[z,x-y]
ISQ /: Times[ISQ[a_List,x_],ISQ[b_List,y_]] := ISQ[a+b,x*y]
ISQ /: Divide[ISQ[a_List,x_],ISQ[b_List,y_]] := ISQ[a-b,x/y]
ISQ /: Power[ISQ[a_List,c_],b_] := ISQ[b*a,Power[c,b]]
ISQ /: Sqrt[ISQ[a_List,c_]] := ISQ[a/2,Sqrt[c]]
ISQ /: Inverse[ISQ[a_List,c_]] := ISQ[-a,1/c]

ISQ /: Times[ISQ[a_List,c_],b_] := ISQ[a,c*b]

(* Measure *)

Measure /: MakeBoxes[Measure[v_, d_ISQ, u_String],StandardForm] := Module[{x=transform[u,d]}, RowBox[{ToBoxes[v], "[", ToBoxes[x], "]",  u}]]

Measure /: Times[Measure[v_,d_,u_],a_] := Measure[v*a,d,u]
Measure /: Plus[Measure[v_,d1,u_],a_] := Measure[v+a,d1,u]
Measure /: Subtract[Measure[v_,d1,u_],a_] := Measure[v-a,d1,u]
Measure /: Subtract[a_,Measure[v_,d1,u_]] := Measure[a-v,d1,u]

Measure /: Times[Measure[v1_,d1_,u_],Measure[v2_,d2_,u_]] := Measure[v1*v2,d1*d2,u]
Measure /: Divide[Measure[v1_,d1_,u_],Measure[v2_,d2_,u_]] := Measure[v1/v2,d1/d2,u]
Measure /: Power[Measure[v_,d_,u_],b_] := Measure[v^b,d^b,u]
Measure /: Plus[Measure[v1_,d1_,u_],Measure[v2_,d2_,u_]] := Measure[v1+v2,d1+d2,u]
Measure /: Subtract[Measure[v1_,d1_,u_],Measure[v2_,d2_,u_]] := Measure[v1-v2,d1-d2,u]
Measure /: Sqrt[Measure[v_,d_,u_]] := Measure[Sqrt[v],Sqrt[d],u]

Measure /: Equal[Measure[v_,ISQ[l_,c_],u_],a_] := And[Total[l]==0,a==v*c]

(*quantity /: Log[b_,Power[B_,_ISQ]] := "test"*)

ratio[d_, u_, s_] := ratio[d,MatchSystem[s,u],UnitSystem[s]]
ratio[d_, u_UnitSystem, s_UnitSystem] := ratiocalc[TRANSFORM[d], u, s]
ratiocalc[ISQ[d_,c_], u_, s_] :=
	BoltzmannConstant[u, s]^d[[1]]*
	ReducedPlanckConstant[u, s]^d[[2]]*
	SpeedOfLight[u, s]^d[[3]]*
	MagneticConstant[u, s]^d[[4]]*
	ElectronMass[u, s]^d[[5]]*
	MonochromaticRadiation540THzLuminousEfficacy[u, s]^d[[6]]*
	MolarMassConstant[u, s]^d[[7]]*
	AngleConstant[u, s]^d[[8]]*
	RationalizationConstant[u, s]^d[[9]]*
	LorentzConstant[u, s]^d[[10]]*
	GravityConstant[u, s]^d[[11]]

Measure /: Entity["UnitSystem", s_][Measure[v_, d_, u_]] := Measure[v*ratio[d, u, s], d, s]
Measure[v_, d_, u_][Entity["UnitSystem",s_]] := Measure[v*ratio[d, u, s], d, s]
Measure[v_, d_, u_][s_String] := Measure[v*ratio[d, u, s], d, s]
ConvertUnit[Measure[v_, d_, u_],s_] := Measure[v*ratio[d, u, s], d, s]

(* 1,2,3,4, 5, 6, 7,  8,9,10,11 *)
(*kB,ħ,𝘤,μ₀,mₑ,Mᵤ,Kcd,A,λ,αL,g₀ *)
(* F,M,L,T, Q, Θ, N,  J,A,Λ, C  *)

TRANSFORM[ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] := ISQ[{-h,l+t+q/2-f-j,3*f+2*h-l-2*t-q/2,q/-2,m+h+n+2*(f+j)-l-t,-n,j,a-2*r,r-q/2,-q-c,l+t-h-2*(f+j)},x]

transform[_String,x_] := transform["Metric",x]

transform["MetricEngineering",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{f,m,l,t,q,h,n,j,a,0,0},x]
transform["GravitationalMetric",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{f+m,0,l-m,t+2*m,q,h,n,j,0,0,0},x]
transform["Metric",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m,f+l,t-2*f,q,h,n,j,0,0,0},x]
transform["Astronomical",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{f,0,l+3*m,t-2*m,q,h,n,j,a,0,0},x]

transform["Gauss",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m+q/2,f+l+(3/2)*q+c,t-2*f-q-c,0,h,0,j,0,0,0},x]
transform["ESU",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m+q/2,f+l+(3/2)*q,t-2*f-q,0,h,0,j,0,0,0},x]
transform["EMU",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m+q/2,f+l+q/2,t-2*f,0,h,0,j,0,0,0},x]
transform["Kennelly",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m+q/2,f+l+q/2,t-2*f,0,h,n,j,0,0,0},x]

transform["Stoney",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m+h,l+t-f,0,q,0,0,0,0,0,0},x]
transform["Electronic",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,0,l+t-f,0,q,0,0,0,0,0,0},x]
transform["QCDoriginal",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,2*f+m+h-l-t,0,0,q,0,0,0,a,0,0},x]
transform["PlanckGauss",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,2*f+m+h-l-t,0,0,0,0,0,0,a,0,0},x]
transform["Planck",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,2*f+m+q+h-l-t,0,0,0,0,0,0,0,0,0},x]
transform["Natural",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,0,0,0,0,0,0,0,0,0,0},x]
transform["NaturalGauss",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,0,0,0,0,0,0,0,a,0,0},x]

transform["Rydberg",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,m+t+h-f,l+2*(t-h)-3*f,0,q,0,0,0,0,0,0},x]
transform["Hartree",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,0,l+2*(t-h)-3*f,0,q,0,0,0,0,0,0},x]
transform["Hubble",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,f+m,l+t-f,0,q,h,n,j,0,0,0},x]
transform["CosmologicalQuantum",ISQ[{f_,m_,l_,t_,q_,h_,n_,j_,a_,r_,c_},x_]] :=
	ISQ[{0,m+2*f-l-t,0,0,q,h,n,j,0,0,0},x]

transform["LorentzHeaviside",x_] := transform["Gauss",x]
transform["Thomson",x_] := transform["EMU",x]
(*transform["Kennelly",x_] := transform["EMU",x]*)
transform["Schrodinger",x_] := transform["Rydberg",x]
transform["QCD",x_] := transform["Planck",x]
transform["QCDGauss",x_] := transform["PlanckGauss",x]
transform["Cosmological",x_] := transform["Hubble",x]

transform["SI2019Engineering",x_] := transform["MetricEngineering",x]
transform["GravitationalSI2019",x_] := transform["GravitationalMetric",x]
transform["British",x_] := transform["GravitationalMetric",x]
transform["British2019",x_] := transform["GravitationalMetric",x]
transform["English",x_] := transform["MetricEngineering",x]
transform["English2019",x_] := transform["MetricEngineering",x]
transform["Survey",x_] := transform["MetricEngineering",x]
transform["Survey2019",x_] := transform["MetricEngineering",x]

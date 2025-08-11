(* ::Package:: *)

masses={1,1};
spins={{0.9,0,0},{-0.9,0,0}};
impactParams={{0,-25,0},{0,25,0}};
\[Beta]={{0.5,0,0},{-0.5,0,0}};


smin=1;
smax=1;
tmin=-1000;tmax=1500;nt=750;
spatialResolution={150,150,150};
size={2000,2000,2000};
origin={-1000,-1000,-1000};
tmintraj=-5000;tmaxtraj=5000;nttraj:=10000;


filepath="Data/"


(* ::Chapter:: *)
(*Parameter definitions*)


metricTensor=SparseArray[{{1,1}->1,{2,2}->-1,{3,3}->-1,{4,4}->-1}];
impactParamsFour=Prepend[#,0]&/@impactParams;
impactParam=impactParamsFour[[2]]-impactParamsFour[[1]];mB=Norm[impactParam];
mB=Norm[impactParam];
lorentzFactors=1/Sqrt[1-(Norm/@\[Beta])^2];
velocities=lorentzFactors*(Prepend[#,1]&/@\[Beta]);
pauliSpinVectors=Table[{lorentzFactors[[bh]]*\[Beta][[bh]] . spins[[bh]],spins[[bh]][[1]]+(lorentzFactors[[bh]]-1)*\[Beta][[bh]][[1]]/Norm[\[Beta][[bh]]]^2 \[Beta][[bh]] . spins[[bh]],spins[[bh]][[2]]+(lorentzFactors[[bh]]-1)*\[Beta][[bh]][[2]]/Norm[\[Beta][[bh]]]^2 \[Beta][[bh]] . spins[[bh]],spins[[bh]][[3]]+(lorentzFactors[[bh]]-1)*\[Beta][[bh]][[3]]/Norm[\[Beta][[bh]]]^2 \[Beta][[bh]] . spins[[bh]]},{bh,2}];
\[Gamma]=lorentzFactors[[1]]*lorentzFactors[[2]]*(1-\[Beta][[1]] . \[Beta][[2]]);(*Lorentz factor of bh 2 if we move into the rest frame of bh 1*)
constantSpinTensor=S[bh_][p1_,p2_]:>TensorContract[TensorProduct[LeviCivitaTensor[4],metricTensor . p1,metricTensor . p2,metricTensor . velocities[[bh]],metricTensor . pauliSpinVectors[[bh]]],{{1,5},{2,6},{3,7},{4,8}}];
replaceParameters={m[bh_]:>masses[[bh]],v[bh_]:>velocities[[bh]],b:>impactParam};
swapHoles={m[bh_]:>m[bh/.{1->2,2->1}],v[bh_]:>v[bh/.{1->2,2->1}],b[bh_]:>b[bh/.{1->2,2->1}],b->-b,\[Tau][bh_]:>\[Tau][bh/.{1->2,2->1}],CE[bh_]:>CE[bh/.{1->2,2->1}],\[Rho]V[bh_]:>\[Rho]V[bh/.{1->2,2->1}],\[Epsilon]V[bh_]:>\[Epsilon]V[bh/.{1->2,2->1}],b\[Epsilon]->-b\[Epsilon],u[bh_]:>u[bh/.{1->2,2->1}],mBt[bh_]:>mBt[bh/.{1->2,2->1}],v[bh_]:>v[bh/.{1->2,2->1}],S[bh_]:>S[bh/.{1->2,2->1}]};
\[Kappa]=Sqrt[32\[Pi]];(*Coupling constant*)


(* ::Chapter:: *)
(*Generate trajectory*)


(* ::Section:: *)
(*Setup*)


gtConstantSpinTensor={S[bh_][p1_]:>TensorContract[TensorProduct[LeviCivitaTensor[4],metricTensor . p1,metricTensor . velocities[[bh]],metricTensor . pauliSpinVectors[[bh]]],{{1,5},{3,6},{4,7}}],constantSpinTensor};
gtReplaceParameters=Join[{\[Tau][bh_]:>\[Tau]/lorentzFactors[[bh]],b[bh_]:>Symbol["cp"<>ToString[bh]]*impactParamsFour[[bh]]+Symbol["cc"<>ToString[bh]]*Prepend[Cross[impactParams[[bh]],\[Beta][[bh]]],0]},replaceParameters];


<<1PM_deflection.wl;
trajectoryBh1Uncorrected[\[Tau]_]=(trajectory/.gtConstantSpinTensor/.gtReplaceParameters)[[2;;4]];
trajectoryBh2Uncorrected[\[Tau]_]=(trajectory/.swapHoles/.gtConstantSpinTensor/.gtReplaceParameters)[[2;;4]];


(* ::Section:: *)
(*Correction*)


(* ::Text:: *)
(*The impact parameter in the trajectory is the impact parameter at \[Tau]=0. Usually the impact parameter is defined as the distance at u->infinity to the impact trajectory. Therefore we have to adjust the impact parameter.*)


tempCorrVar={cp1->tcp1,cc1->tcc1,cp2->tcp2,cc2->tcc2};
limits=Join[{\[Tau]->-Infinity},tempCorrVar];
correction=FindRoot[{Limit[trajectoryBh1Uncorrected[\[Tau]] . impactParams[[1]],limits]==(Norm[impactParams[[1]]]^2),Limit[trajectoryBh1Uncorrected[\[Tau]] . Cross[impactParams[[1]],\[Beta][[1]]],limits]==0,Limit[trajectoryBh2Uncorrected[\[Tau]] . impactParams[[2]],limits]==Norm[impactParams[[2]]]^2,Limit[trajectoryBh2Uncorrected[\[Tau]] . Cross[impactParams[[2]],\[Beta][[2]]],limits]==0},{{tcp1,1},{tcc1,0},{tcp2,1},{tcc2,0}}];
trajectoryBh1[\[Tau]_]=trajectoryBh1Uncorrected[\[Tau]]/.tempCorrVar/.correction;
trajectoryBh1Comp=With[{code=N[trajectoryBh1[\[Tau]]]},Compile[{{\[Tau],_Real}},code,CompilationTarget->"C",Parallelization->True,RuntimeAttributes->Listable]];
trajectoryBh2[\[Tau]_]=trajectoryBh2Uncorrected[\[Tau]]/.tempCorrVar/.correction;
trajectoryBh2Comp=With[{code=N[trajectoryBh2[\[Tau]]]},Compile[{{\[Tau],_Real}},code,CompilationTarget->"C",Parallelization->True,RuntimeAttributes->Listable]];


(* ::Section:: *)
(*Make ready for export*)


urangetraj=Range[tmintraj,tmaxtraj,(tmaxtraj-tmintraj)/(nttraj-1)];
trajectoryBh1Data=Join[ArrayReshape[urangetraj,{nttraj,1}],trajectoryBh1Comp[urangetraj],2];
spinBh1Data=Join[ArrayReshape[urangetraj,{nttraj,1}],ConstantArray[spins[[1]],Length[urangetraj]],2];
trajectoryBh2Data=Join[ArrayReshape[urangetraj,{nttraj,1}],trajectoryBh2Comp[urangetraj],2];
spinBh2Data=Join[ArrayReshape[urangetraj,{nttraj,1}],ConstantArray[spins[[2]],Length[urangetraj]],2];
trajectoryData={"/trajectory1/position.dat"->{"Data"->trajectoryBh1Data,"DataFormat"->"Real32"},"/trajectory2/position.dat"->{"Data"->trajectoryBh2Data,"DataFormat"->"Real32"},"/trajectory1/spin.dat"->{"Data"->spinBh1Data,"DataFormat"->"Real32"},"/trajectory2/spin.dat"->{"Data"->spinBh2Data,"DataFormat"->"Real32"}};


(* ::Chapter:: *)
(*Generate waveform*)


(* ::Section:: *)
(*Setup*)


dot[p1_,p2_]:=TensorContract[TensorProduct[SparseArray[{{1,1}->1,{2,2}->-1,{3,3}->-1,{4,4}->-1}],p1,p2],{{1,3},{2,4}}];
P=masses[[1]]*velocities[[1]]+masses[[2]]*velocities[[2]];
tE=Sqrt[dot[P,P]];
V=P/tE;
VO=(velocities[[1]]*dot[velocities[[2]],V]-velocities[[2]]*dot[velocities[[1]],V])/Sqrt[\[Gamma]^2-1];
L=(masses[[1]]*masses[[2]])/tE TensorContract[TensorProduct[LeviCivitaTensor[4],impactParam,velocities[[1]],velocities[[2]]],{{2,5},{3,6},{4,7}}];
nb=impactParam/Sqrt[-dot[impactParam,impactParam]];
nL=L/Sqrt[-dot[L,L]];
\[Rho]=V+VO*Cos[\[Theta]]+(nb*Cos[\[Phi]]+nL*Sin[\[Phi]])*Sin[\[Theta]];
unit\[Theta]=-VO*Sin[\[Theta]]+(nb*Cos[\[Phi]]+nL*Sin[\[Phi]])*Cos[\[Theta]];
unit\[Phi]=-nb*Sin[\[Phi]]+nL*Cos[\[Phi]];
\[Epsilon]=unit\[Theta]+I*unit\[Phi];
\[Rho]timesv=Table[dot[\[Rho],velocities[[bh]]],{bh,2}]
\[Epsilon]timesv=Table[dot[\[Epsilon],velocities[[bh]]],{bh,2}];
btimes\[Epsilon]=dot[impactParam,\[Epsilon]];
retardedTime=Table[(u-dot[\[Rho],impactParamsFour[[bh]]])/\[Rho]timesv[[bh]],{bh,2}];
shiftetImpactCOM=Sqrt[mB^2-retardedTime[[1]]^2-retardedTime[[2]]^2+2\[Gamma]*retardedTime[[1]]*retardedTime[[2]]];
shiftetImpactRest=Table[Sqrt[mB^2+(\[Gamma]^2-1)retardedTime[[bh/.{1->2,2->1}]]^2],{bh,2}];
gwReplaceParameters=Join[{\[Rho]V[bh_]:>\[Rho]timesv[[bh]],\[Epsilon]V[bh_]:>\[Epsilon]timesv[[bh]],b\[Epsilon]->btimes\[Epsilon],u[bh_]:>retardedTime[[bh]],mBt[0]->shiftetImpactCOM,mBt[bh_]:>shiftetImpactRest[[bh]]},replaceParameters];


<<ancillary.m
halfWaveform = m[1]m[2]Sum[(Symbol["\[Alpha]1"<>ToString[s]]+Divide[Symbol["\[Beta]1"<>ToString[s]],mBt[0]^(2s+2)])/mBt[1]^(2s+1),{s,smin,smax}];
waveform[u_,\[Theta]_,\[Phi]_]=1/\[Kappa] (halfWaveform +(halfWaveform/.swapHoles)/.constantSpinTensor/.gwReplaceParameters);
waveformComp=With[{code=N[waveform[u,\[Theta],\[Phi]]]},Compile[{{u,_Real},{\[Theta],_Real},{\[Phi],_Real}},code,CompilationTarget->"C",Parallelization->True,RuntimeAttributes->Listable]];


(* ::Section:: *)
(*Generate waveform*)


Nspatial = spatialResolution[[1]] spatialResolution[[2]] spatialResolution[[3]];
trange = N[Range[tmin, tmax, (tmax - tmin)/(nt - 1)]];
step := size/(spatialResolution - 1);
end := origin + size;
cartesianCoords := Flatten[N[Table[{k, j, i}, {i, origin[[1]], end[[1]], step[[1]]}, {j, origin[[2]], end[[2]], step[[2]]}, {k, origin[[2]], end[[2]], step[[2]]}]], 2];
coordinates = N[ToSphericalCoordinates[cartesianCoords]];


(* ::Section:: *)
(*Make ready for export*)


timestepData={"/waveform/timesteps.dat"->{"Data"->trange,"DataFormat"->"Real32"}};


(* ::Chapter:: *)
(*Export as HDF file*)


Export[filepath <> "trajectory1.csv", N[trajectoryBh1Data], "CSV"];
Export[filepath <> "trajectory2.csv", N[trajectoryBh2Data], "CSV"];
hdfFile = filepath <> "gravitationalBremsstrahlung.h5";
Export[hdfFile,{"Datasets"->Join[trajectoryData,timestepData],"Attributes"->{"/"->{"m[1]"->masses[[1]],"m[2]"->masses[[2]],"\[Beta][1]"->\[Beta][[1]],"\[Beta][2]"->\[Beta][[2]],"b[1]"->impactParams[[1]],"b[2]"->impactParams[[2]],"s[1]"->spins[[1]],"s[2]"->spins[[2]],"smin"->smin,"smax"->smax},"/waveform"->{"tmin"->tmin,"tmax"->tmax,"nt"->nt,"dim_x"->spatialResolution[[1]],"dim_y"->spatialResolution[[2]],"dim_z"->spatialResolution[[3]],"origin_x"->origin[[1]],"origin_y"->origin[[2]],"origin_z"->origin[[3]],"size_x"->size[[1]],"size_y"->size[[2]],"size_z"->size[[3]]}}},"Rules"];
For[i = 1, i <= nt, i++, utable := N[ConstantArray[trange[[i]], Nspatial] - coordinates[[All, 1]]]; waveformGridOverTime := waveformComp[utable, coordinates[[All, 2]], coordinates[[All, 3]]]/coordinates[[All, 1]]; Export[hdfFile,{"/waveform/strain/"<>ToString[i]<>".dat"->{"Data"->ReIm[waveformGridOverTime],"DataFormat"->"Real32"}},OverwriteTarget->"Append"];Run[ "echo -n '" <> ToString[i]<>"/"<> ToString[nt] <> "\r'"]]

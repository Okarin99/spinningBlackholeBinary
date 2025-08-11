(* ::Package:: *)

(* ::Title:: *)
(*1PM trajectory of two spinning black holes*)


(* ::Section:: *)
(*Notation*)


(* ::Text:: *)
(*-- b(relative impact parameter b^mu=b_2^mu-b_1^mu (at \[Tau][1] =0))*)
(*-- mB (modulus of the impact parameter (at \[Tau][1] =0))*)
(*-- S[i][\[Mu],\[Nu]] := S_i^{\[Mu]\[Nu]} (constant spin tensor)*)
(*-- v[i] (initial velocities of BHs)*)
(*-- \[Kappa] (gravitational coupling \[Kappa]^2 = 32 pi G)*)
(*-- m[i] (masses of the BHs)*)
(*-- \[Gamma] (relative lorentz boost factor of BH 2 in the frame of BH 1)*)
(*-- \[Tau][1] (proper time in the frame of BH 1)*)


(* ::Section:: *)
(*Tajectory*)


(* ::Subsection:: *)
(*z^mu(\[Tau][1]) in zero order in spin*)


defl0spin = -(1/(32 mB^2 \[Pi] (-1+\[Gamma]^2)^(3/2)))m[2]\[Kappa]^2 ((-1+2 \[Gamma]^2) (\[Tau][1]-\[Gamma]^2 \[Tau][1]+Sqrt[-1+\[Gamma]^2] (mB-Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2])) b -
mB^2 Log[(Sqrt[-1+\[Gamma]^2] \[Tau][1]+Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2])/mB] (v[1]+\[Gamma] (-3+2 \[Gamma]^2) v[2]));


(* ::Subsection:: *)
(*z^mu(\[Tau][1]) in linear order in spin*)


defllspin = (m[2] \[Kappa]^2 (-2 mB \[Gamma] (mB^2 Sqrt[-1+\[Gamma]^2]-mB Sqrt[-1+\[Gamma]^2] Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2]+2 (-1+\[Gamma]^2) \[Tau][1] (Sqrt[-1+\[Gamma]^2] \[Tau][1]+Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2])) b (S[1][b,v[2]]-S[2][b,v[1]]) +
mB^3 ((-1+\[Gamma]^2)^(3/2) \[Tau][1] S[1][b]+\[Gamma] (mB^2 Sqrt[-1+\[Gamma]^2]-mB Sqrt[-1+\[Gamma]^2] Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2]+2 (-1+\[Gamma]^2) \[Tau][1] (Sqrt[-1+\[Gamma]^2] \[Tau][1]+Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2])) S[1][v[2]] -
2 \[Tau][1] (-\[Gamma] (-1+\[Gamma]^2)^(3/2) S[2][b]+\[Gamma] (-1+\[Gamma]^2) (Sqrt[-1+\[Gamma]^2] \[Tau][1]+Sqrt[mB^2+(-1+\[Gamma]^2) \[Tau][1]^2]) S[2][v[1]]+Sqrt[-1+\[Gamma]^2] (S[1][b,v[2]]-S[2][b,v[1]]) (\[Gamma] v[1]-v[2]))))) /
(32 mB^4 \[Pi] (-1+\[Gamma]^2)^(3/2) Sqrt[mB^4+mB^2 (-1+\[Gamma]^2) \[Tau][1]^2]);


(* ::Subsection:: *)
(*Trajectory*)


trajectory = b[1]+v[1]\[Tau][1]+defl0spin+defllspin

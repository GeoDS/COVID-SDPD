 function results =lr_f_err_c(y,x,W,W2,N)
%PURPOSE : Computes the LR statistic to test the presence of an endogenous spatial lag 
%when spatially autocorrelated errors are already accounted for in a fixed
%effects panel data model (with individual effects only, no time effects).
%Under the null, lamdba can be different from zero while rho is set to
%zero. Under the alternative, both spatial parameters are not constrained

% y : (NT x 1) vector of the dependent variable
% x : (NT x K) matrix of independent variables (WITHOUT CONSTANT)
% W : (N x N) weight matrix for the autoregressive part  
% W2 : (N x N) weight matrix for spatially autoregressive errors 
% Note :  the code has been written assuming possibly different Weight matrices for both the
% SAR and SEM parts.
% N : number of individuals.

% Data organisation : i is fast moving index and t is slow moving index (all 
% individuals  for period 1, all individuals for period 2, ...). 
% If it is not the case,apply the function trans_tslow.m 

% The input of the test must be the original data. The Lee and Yu
% transformation is applied inside the test.

% Code has been programmed for balanced panel
%--------------------
% Written by N. Debarsy* and C. Ertur**
% * University of Namur
%   Centre de recherches en Economie R�gionale et Politique Economique (CERPE)
%   Rempart de la vierge, 8
%   5000 Namur, Belgium
%   nicolas.debarsy@fundp.ac.be

%** Universit� d'Orl�ans
%   UFR Droit-Economie-Gestion
%   Laboratoire d'Economie d'Orl�ans - UMR 6221 CNRS
%   Domaine Universitaire
%   Rue de Blois - BP 6739
%   45067 ORLEANS Cedex 2, France
%   cem.ertur@univ-orleans.fr

% REFERENCE :
% Debarsy N, C. Ertur (2009), Testing for spatial autocorrelation in a fixed effects
% panel data model, Document de recherche du LEO, 2009-12. 
% --------------------------------------------------------

[nt k]=size(x);
% Computation of the unconstrained model (SAC estimation)
info.lflag=0;
res1 = sarar_panel_FE_LY(y,x,W,W2,N,info);
liku=res1.lik;

% Computation of the constrained model (SEM model) 
res=sem_panel_FE_LY(y,x,W2,N,info);
likr=res.lik;

% Computation of the test
lr=2*(liku-likr);
prob=1-chis_prb(lr,1);

results.meth='lr_f_err_c';
results.lr=lr;
results.prob=prob;
results.Wsar= inputname(3);
results.Wsem=inputname(4);
results.nvar = k;
results.chi_1= 6.635;
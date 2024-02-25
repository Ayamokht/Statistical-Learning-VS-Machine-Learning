
/***********************DGP avec multicolinéarité****************/


proc iml;
call randseed(1234);

N = 500;
nb_var_noncorr = 45;
nb_var_noncorr= 5;
cov_noncorr=I(45);
mean_noncorr = j(45,1,0);
it=1000;
/*varl = ("Intercept"||('X1':'X50'))`;
indpro=j(nrow(varl),1,0);*/ /*pour les graphiques avec les probabilités individuelles*/


goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*décomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'};   /*vraies variables utilisées dans le DGP*/


/** 1 AIC **/

/****************** PARTIE 1 ******************/

do a=1 to it;

/** 45 variables non-corélées **/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/*coefficients de correlation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);


X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);


/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;



proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=AIC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;


endsubmit;  

use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;


end;
colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_AIC=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_AIC[colname=colres];

create table_AIC from res_AIC[colname=colres];
append from res_AIC;
close table_AIC;



/** 2 AICC **/


/****************** PARTIE 1  ******************/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisées dans le DGP*/

do a=1 to it;

/** 45 variables non-corrélées **/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/*coefficients de corrélation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);

X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);

/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;


proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=AICC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  

use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);



/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;

end;
colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_AICC=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_AICC[colname=colres];

create table_AICC from res_AICC[colname=colres];
append from res_AICC;
close table_AICC;




/** 3 BIC **/

/****************** PARTIE 1 ******************/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*décomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisées dans le DGP*/


do a=1 to it;

/** 45 variables non-corélées**/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/*coefficients de corrélation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);


X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);

/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;


proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=BIC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  


use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);



/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;


end;
colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_BIC=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_BIC[colname=colres];

create table_BIC from res_BIC[colname=colres];
append from res_BIC;
close table_BIC;



/*** 4 SBC **/

/****************** PARTIE 1 ******************/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*décomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisées dans le DGP*/


do a=1 to it;

/** 45 variables non-corrélées **/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/*coefficients de corrélation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);


X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);

/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;


proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=SBC stop=PRESS lscoeffs);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  


use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);



/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;



end;

colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_SBC=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_SBC[colname=colres];

create table_SBC from res_SBC[colname=colres];
append from res_SBC;
close table_SBC;


/** 5 CP **/

/****************** PARTIE 1 ******************/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*décomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisées dans le DGP*/


do a=1 to it;

/** 45 variables non-corélées **/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/* coefficients de corrélation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);


X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);


/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;



proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=CP stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  


use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);



/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;


end;
colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_CP=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_CP[colname=colres];

create table_CP from res_CP[colname=colres];
append from res_CP;
close table_CP;


/** 6 CV **/

/****************** PARTIE 1 ******************/


goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*décomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisées dans le DGP*/

do a=1 to it;

/** 45 variables non-corrélées **/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/* coefficients de corrélation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);


X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);

/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;


proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=CV stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  


use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);



/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;


end;
colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_CV=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_CV[colname=colres];

create table_CV from res_CV[colname=colres];
append from res_CV;
close table_CV;



/**7 PRESS**/

/****************** PARTIE 1 ******************/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*décomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisées dans le DGP*/

do a=1 to it;

/** 45 variables non-corrélées **/

x_noncorr = randnormal(N, mean_noncorr, I(45));

/* coefficients de corrélation */
coeff= {1 0.8 0.6 0.4 0.5};
mat_T= toeplitz(coeff);


X_corr=randnormal(N,{0,0,0,0,0},mat_T);
c=corr(X_corr);


/* concaténation des deux */
X= X_corr || X_noncorr; 

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(500,1,2))*0.05;

resu1=Y||X;
varnames="y"||("X1":"X50"); 

create colinearity_table from resu1[colnames=varnames];
append from resu1;
close colinearity_table;

/****************** PARTIE 2 ******************/

submit;


proc glmselect data= colinearity_table noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=PRESS stop=SBC LSCOEFFS);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  



use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enlève y */
varnames= ncol(Mat);



/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnées
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionnées*/

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail général*/
	

/*II. Décomposition de fail en plusieurs sous-cas*/

/*2.1 almost_gf : il retourne 5 variables mais pas toutes sont bonnes (en dehors de la cst), ex : intercept,X1,X2,X3,X4+X45, cardmat=6, cardintersection=5<6 */
if cardmat=6 & cardintersection<6 then
		almost_gf=almost_gf+1;
	else almost_gf=almost_gf;


/*2.2 overfit_a : il retourne >5 variables mais TOUTES les 5 bonnes sont dedans, ex : intercept,X1,X2,X3,X4,X5+X9,X10, cardmat=8>6, cardintersection=6*/
if cardmat>6 & cardintersection=6 then 
		overfit_a=overfit_a+1;
	else overfit_a=overfit_a;

/*2.3 overfit_b : il retourne plus que 5 variables, mais PAS TOUTES les bonnes sont dedans, ex : intercept,X1,X2,X3+X7,X8,X9, cardmat=7>6, cardintersection=4<6*/
if cardmat>6 & cardintersection<6 then
		overfit_b=overfit_b+1;
	else overfit_b=overfit_b;

	/*2.4 underfit_a : il retourne moins que 5 variables, mais toutes bonnes, ex : intercept,X1,X2,X3, cardmat=4<6, cardintersection=cardmat=4*/
if cardmat<6 & cardmat=cardintersection then
		underfit_a=underfit_a+1;
	else underfit_a=underfit_a;

 	/*2.5 underfit_b : il retourne moins que 5 variables, pas toutes sont bonnes, ex : intercept,X1,X2+X7,X8, cardmat=5, cardintersection=3<5*/
if cardmat<6 & cardintersection<cardmat then
		underfit_b=underfit_b+1;
	else underfit_b=underfit_b;


/*Pour les probabilités individuelles
do i=1 to nrow(Mat);
l=loc(varl=Mat[i]);
if ncol(l)^=0 then indpro[l]=indpro[l]+1;
else indpro[l]=indpro[l];


end;
end;
indpro=indpro/(a-1);
print indpro[rowname=varl label={'Individual probas'}];*/

end;
colres = {"goodfit", "fail", "almost_gf", "overfit_a", "overfit_b", "underfit_a", "underfit_b"};
res_PRESS=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_PRESS[colname=colres];

create table_PRESS from res_PRESS[colname=colres];
append from res_PRESS;
close table_PRESS;


data combined_dataset;

length Criterion $20; 
Criterion = "AIC";
it=1000;
 
   set table_AIC table_AICC table_BIC table_SBC table_CP table_CV table_PRESS;

   if _n_ = 2 then Criterion = "AICC";
   if _n_ = 3 then Criterion = "BIC";
   if _n_ = 4 then Criterion = "SBC";
   if _n_ = 5 then Criterion = "CP";
   if _n_ = 6 then Criterion = "CV";
   if _n_ = 7 then Criterion = "PRESS";


run;


data combined_dataset ; set combined_dataset;
p_good_fit = (goodfit / it) * 100;
p_almost_good_fit= (almost_gf / it) * 100;
p_good_overfit = (overfit_a / it) * 100;
p_bad_overfit = (overfit_b / it) * 100;
p_good_underfit = (underfit_a / it) * 100;
p_bad_underfit = (underfit_b / it) * 100;
run;
quit; 

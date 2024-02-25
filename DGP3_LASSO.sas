
/****************************** DGP3 avec valeurs extremes  ****************************/



/*1. AIC */
proc iml;
call randseed(1234);
Mean=j(1,50,0);
Cov=I(50);
N=500;
boucle=1000;

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il ne rate pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'}; /*vraies variables utilisees dans le DGP*/

do a=1 to boucle;
X=randnormal(N, Mean, Cov);
/*on genere un vecteur avec des valeurs entre 0 et 1*/


do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; /* Cr�ation Table SAS pour DGP_outliers */
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;

/* crit�re statistique: forward backword stepwise */

proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=AIC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  


/**** turn the table DGP1 into matrix *****/
use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* On enl�ve Y */
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, nombre de variables communes */
cardmat=nrow(Mat); /*card = nombre de variables qui ont �t� selectionn�es */

/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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



/* 2. AICC*/

/*************** PARTIE 1 **************/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'};


do a=1 to boucle;
X=randnormal(N, Mean, Cov);

do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; 
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;

proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=AICC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;


endsubmit;  



use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1];
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, combien de variables communes ?*/
cardmat=nrow(Mat); /*card = combien de variables ont ete selectionnes ?*/


/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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



/*3. BIC*/


goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'};

do a=1 to boucle;
X=randnormal(N, Mean, Cov);



do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; 
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;



proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=BIC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  



use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* drop the column y */
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, combien de variables communes ?*/
cardmat=nrow(Mat); /*card = combien de variables ont ete selectionnes ?*/


/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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


/*4. SBC*/


goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/
true={'Intercept','X1','X2','X3','X4','X5'};


do a=1 to boucle;
X=randnormal(N, Mean, Cov);


do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; 
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;



proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=SBC stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;

endsubmit;  


use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; 
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionn�es*/


/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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



/*5. CP*/

goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'};

do a=1 to boucle;
X=randnormal(N, Mean, Cov);



do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; 
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;



proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=CP stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  


use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; 
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, combien de variables communes ?*/
cardmat=nrow(Mat); /*card = combien de variables ont ete selectionnes ?*/


/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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


/*6. CV*/


goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'};

do a=1 to boucle;
X=randnormal(N, Mean, Cov);



do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; 
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;



proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=CV stop=SBC);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  



use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enl�ve Y */
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionn�es*/


/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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


/*7. PRESS*/
goodfit=0; /*il retourne les 5 bonnes variables*/
fail=0; /*il rate=pas goodfit*/

/*decomposition du fail*/

underfit_a=0; /*il retourne moins que 5 variables, toutes bonnes*/
underfit_b=0;/*il retourne moins que 5 variables, PAS toutes bonnes*/
overfit_a=0; /*il retourne plus que 5 variables, toutes les bonnes sont dedans*/
overfit_b=0;/*il retourne plus que 5 variables, PAS toutes les bonnes sont dedans*/
almost_gf=0; /*il retourne 5 mais une combinaison de bonnes et de mauvaises*/

true={'Intercept','X1','X2','X3','X4','X5'};

do a=1 to boucle;
X=randnormal(N, Mean, Cov);



do i = 1 to N;
do t=1 to 3;
u=uniform(0);
    if u > 0.95 then
        X[i, t] = rand('Normal', 4, 1);
    else
        X[i, t] = rand('Normal', 0, 1);
end;
end;

Y=0.8*X[,1]-0.7*X[,2]+0.5*X[,3]-1.2*X[,4]+0.9*X[,5]+normal(j(N,1,2))*0.05;

varnames="y"||("X1":"X50"); 
resu=Y||X;
create table_outliers from resu[colnames=varNames]; 
append from resu;
close table_outliers;


/****************** PARTIE 2 ******************/

submit;



proc glmselect data= table_outliers noprint outdesign= DGP2_results;
model Y =X1-X50 / selection=lasso(choose=PRESS stop=SBC LSCOEFFS);
run;
proc transpose data=DGP2_results out=DGP2_results; run;



endsubmit;  



use DGP2_results;
read all;
close DGP2_results;
Mat=_Label_[1:nrow(_Label_)-1]; /* on enl�ve Y */
varnames= ncol(Mat);


/*********** PARTIE 3: LES METRIQUES ***************/


/*mat-les variables selectionnees
true-les vraies variables X1-X5*/


intersection=xsect(true,Mat);
cardintersection=ncol(intersection); /*card de l'intersection, le nombre de variables communes */
cardmat=nrow(Mat); /*card = le nombre de variables selectionn�es*/


/*I.goodfit vs fail*/
if cardmat=6 & cardintersection=6 then
		goodfit=goodfit+1;
	else fail=fail+1; /*fail general*/
	

/*II. Decomposition de fail en plusieurs sous-cas*/

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
res_PRESS=goodfit||fail||almost_gf||overfit_a||overfit_b||underfit_a||underfit_b;
print res_PRESS[colname=colres];

create table_PRESS from res_PRESS[colname=colres];
append from res_PRESS;
close table_PRESS;



data combined_dataset;

length Criterion $20; 
Criterion = "AIC";
boucle=1000;
 
   set table_AIC table_AICC table_BIC table_SBC table_CP table_CV table_PRESS;

   if _n_ = 2 then Criterion = "AICC";
   if _n_ = 3 then Criterion = "BIC";
   if _n_ = 4 then Criterion = "SBC";
   if _n_ = 5 then Criterion = "CP";
   if _n_ = 6 then Criterion = "CV";
   if _n_ = 7 then Criterion = "PRESS";

run;


data combined_dataset ; set combined_dataset;
p_good_fit = (goodfit / boucle) * 100;
p_almost_good_fit= (almost_gf / boucle) * 100;
p_good_overfit = (overfit_a / boucle) * 100;
p_bad_overfit = (overfit_b / boucle) * 100;
p_good_underfit = (underfit_a / boucle) * 100;
p_bad_underfit = (underfit_b / boucle) * 100;
run;

quit;


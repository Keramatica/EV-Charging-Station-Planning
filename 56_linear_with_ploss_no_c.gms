
************** Linearization Approach *****************

set
e   Index of piecewises /e1*e10/;

*****************************

Scalars
ee    The number of piecewises /10/,
ilim_max /0.1/,
dilim_max  /0.01/,
ilr_max /0.1/,
dilr_max  /0.01/;
******************************
          set i index for bus of DN /1*6/  ;
         Alias(i,j) ;
         set  k index for harmonics/1*11/
sets
         z index for  location of charging station /1*5/
         slack(i) index for couped bus -upstream /1/
         nl index for nonlinear load /sp/
         a index for level of operation/1*4/
h index for harmonic orders/1,5,7,11,13,17,19,23,25,29,31/
*h index for harmonic orders/1,5/

         r index for o-d pair /1-3,1-4,2-3,2-4/
         LLoad(i)/2,3,4,5,6/
         node_bus(z,i)  cuopling  Transport and Dist System /1.3,2.6,3.4,4.5,5.2/
         nl_bus(i,nl)   nonlinear load and corresponding bus/5.sp/
         line(i,j)  Dist_line/1.2,2.3,2.4,4.5,5.6/
         s index for season of year /B-P,T,Z/
         ti index for time of day /1*24/
       d1/1*16/
d2/1*8/
d4/1*4/ ;

*^^^^^^^^^^^^^^^^^^^^^^^
sets
        an link in transportation network/a1*a10/
        z_an(z,r,an) /1.1-3.a2,1.1-3.a6,2.1-3.a2,2.1-3.a4,2.1-3.a9,2.1-3.a6,3.1-3.a2,3.1-3.a6,4.1-3.a2,4.1-3.a10,4.1-3.a8,5.1-3.a2,5.1-3.a6
,1.1-4.a2,1.1-4.a10,2.1-4.a2,2.1-4.a4,2.1-4.a9,2.1-4.a10,3.1-4.a2,3.1-4.a6,3.1-4.a3,4.1-4.a2,4.1-4.a10,5.1-4.a2,5.1-4.a10,1.2-3.a9,1.2-3.a7,1.2-3.a1
,2.2-3.a9,2.2-3.a6,3.2-3.a9,3.2-3.a6
,4.2-3.a9,4.2-3.a10,4.2-3.a8,5.2-3.a9,5.2-3.a6,1.2-4.a9,1.2-4.a7,1.2-4.a2,1.2-4.a10,2.2-4.a9,2.2-4.a10,3.2-4.a9,3.2-4.a6,3.2-4.a3,4.2-4.a9,4.2-4.a10
,5.2-4.a9,5.2-4.a10/
;
parameter
ca(an) capacity of each link/a1 23,a2 21,a3 22,a4 14,a5 30,a6 20,a7 21,a8 18,a9 24,a10 17/
FFT(an) free flow travel time/a1 15,a2 8,a3 14,a4 7,a5 20,a6 5,a7 11,a8 11,a9 10,a10 9/
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ;


parameters
         c_fix(z) fix costof location cs ($) /1 163000,2  163000,3  163000,4  163000,5  163000/
         cs_p(z) cost for one charge point /1 20907.44,2 20907.44,3 20907.44,4 20907.44,5 20907.44/
         psp(i) active power for one charge point /1 60,2 60,3 60,4 60,5 60,6 60/
         cs(i) cost of ac-dc converter $.pu(100mva)/1 63220,2 63220,3 63220,4 63220,5 63220,6 63220/

pg(i) generation of active power from power flow /1 0,2 0,3 0,4 0,5 0,6 0/
qg(i) generation of reactive power from power flow /1 0,2 0.015,3 0,4 0,5 0,6 0/
 pdo(i) active power for nomrmal load /1 0,2 0.00522,3 0,4 0.00882,5 0,6 0.00882/
 pnl(nl) active power for nomrmal load   /sp 0.00527/
qnl(nl) active power for nomrmal load /sp 0.00242/
qdo(i) active power for nomrmal load  /1 0,2 0.00174,3 0,4 0.00994,5 0,6 0.00994/

cp  cost of active power duringg 24 hours $.pu(mva)  /125/
cq  cost of active power duringg 24 hours $.pu(mva)  /18.75/

D(r) /1-3 210,1-4 185,2-3 150,2-4 130/
k1(d1)/1 1,2 2,3 3,4 4,5 5,6 6,7 7,8 8,9 9,10 10,11 11,12 12,13 13,14 14,15 15,16 16/
k2(d2)/1 1,2 2,3 3,4 4,5 5,6 6,7 7,8 8/
k4(d4)/1 1,2 2,3 3,4 4/
n/16/

         scalar bzde /0.9/
         scalar beta coefient for time travel to cost ($.h)/8.22/
         scalar cu cunsumption of PEV (KWHdivKM) /0.15/
         scalar yhat maximum number of charge point /4/
         scalar deltat between 2 time /1/
         scalar rng driving range km /100/
         scalar THDmax /0.05/
         scalar IHDmax /0.03/
         SCALAR ILMAX/1.05/
         scalar Sbase MVA /100/
         scalar vbase KV/23/

         scalar vrmsmin lowerest voltage /0.95/
         scalar vrmsmax   highest volrage /1.05/

M/100000/ ;
set ii /1*26/;
alias(ii,jj);

table Voltage(ii,*)'Voltage Mag'
              Real         Image

*** answered***
1              -0.18        -1
2              -0.175      -0.99
3              -0.17       -0.98
4              -0.165      -0.97
5              -0.16       -0.96
6              -0.155      -0.95
7              -0.15       -0.94
8              -0.145      -0.93
9              -0.14       -0.92
10             -0.135      -0.91
11             -0.13       -0.9
12             -0.125      -0.88
13             -0.12       -0.87
14             -0.115      -0.86
15             -0.11       -0.85
16             -0.1        -0.84
17             -0.09       -0.83
18             -0.089      -0.82
19             -0.088      -0.81
20             -0.087      -0.8
21             -0.086      -0.79
22             -0.085      -0.78
23             -0.084      -0.77
24             -0.083      -0.76
25             -0.082      -0.75
26             -0.081      -0.74

;

table Current(jj,*)'Current Mag'

               Real                      Image

***answered**          1        -0.960000000000000        -0.960000000000000
1        -0.190000000000000        -0.0900000000000000
2        -0.180000000000000        -0.0800000000000000
3        -0.170000000000000        -0.0700000000000000
4        -0.160000000000000        -0.0600000000000000
5        -0.150000000000000        -0.0500000000000000
6        -0.140000000000000        -0.0400000000000000
7        -0.130000000000000        -0.0300000000000000
8        -0.120000000000000        -0.0200000000000000
9        -0.110000000000000        -0.0100000000000000
10       -0.100000000000000        0
11       -0.090000000000000        0.0100000000000000
12       -0.080000000000000        0.0200000000000000
13       -0.070000000000000        0.0300000000000000
14       -0.0600000000000000       0.0400000000000000
15       -0.0500000000000000       0.0500000000000000
16       -0.0400000000000000       0.0600000000000000
17       -0.0300000000000000       0.0700000000000000
18       -0.0200000000000000       0.0800000000000000
19       -0.0100000000000000       0.0900000000000000
20        0                         0.100000000000000
21        0.00999999999999998       0.110000000000000
22        0.0200000000000000        0.120000000000000
23        0.0300000000000000        0.130000000000000
24        0.0400000000000000        0.140000000000000
25        0.0500000000000000        0.150000000000000
26        0.0600000000000000        0.160000000000000 ;

parameters
ratio
x_h
bb_h
b_h
g_h
y_h;
                                                                                                                                                                                                                                        ;


table chr(nl,k) harmonic current factor
         1       2       3       4       5       6       7       8       9      10      11
sp       1     0.19    0.131   0.072   0.056   0.033   0.024   0.012   0.008   0.002   0.002;
table chim(nl,k) harmonic current factor
         1       2       3       4       5       6       7       8       9      10      11
sp       1     0.19    0.131   0.072   0.056   0.033   0.024   0.012   0.008   0.002   0.002;




Table   branch(i,j,*)  in pu

                 r            x            b       rateA    rateB   rateC     tap      an        st        min         max


1.2            0.0004        0.0069        0        0        0        0        0        0        1        -360        360
2.3            0.0527        0.0028        0        0        0        0        0        0        1        -360        360
2.4            0.0527        0.0028        0        0        0        0        0        0        1        -360        360
4.5            0.0527        0.0028        0        0        0        0        0        0        1        -360        360
5.6            0.2595        0.1462        0        0        0        0        0        0        1        -360        360 ;



table deltan(z,a)
         1       2       3       4
1        1       1       1       1
2        1       1       1       1
3        1       1       1       1
4        1       1       1       1
5        1       1       1       1 ;
table c_wait(z,a)  cost of waiting because of state of level of opreation charging station  ($.h)
         1       2       3       4
1        1       1.5     3       6
2        1       1.5     3       6
3        1       1.5     3       6
4        1       1.5     3       6
5        1       1.5     3       6;


table TR(r,z)
        1       2       3      4       5
1-3     13      30      13     28      13
1-4     17      34      27     17      17
2-3     34      15      15     30      15
2-4     38      19      29     19      19 ;
table co_p(s,ti) coeffiecent for active power
       1      2       3       4       5       6       7       8       9       10      11      12      13      14     15       16      17      18        19      20      21      22      23      24
B-P    0.45   0.45    0.45    0.45    0.45    0.45    0.45    0.55    0.55    0.55    0.55    0.6     0.6     0.6    0.6      0.6     0.6     0.6      0.6     0.6     0.6     0.6     0.6     0.6
T      0.6    0.6     0.6     0.6     0.6     0.6     0.6     0.8     0.8     0.8     0.8     0.8     1       1       1        1       1      0.9      0.9     0.9     0.9     0.8     0.8     0.8
Z      0.7    0.7     0.7     0.7     0.7     0.7     0.7     0.82    0.82    0.82    0.82    0.82    0.73    0.73   0.73     0.7     0.73    0.73     0.73    0.8     0.8     0.75    0.75    0.75
;

table co_price(s,ti) coeffiecent for price of electricity
       1    2      3       4       5       6       7       8       9        10      11      12      13      14      15      16      17      18       19      20      21      22      23      24
B-P    40  40      40      40      40      60      60      60      60       60      50      50      50      50      50      40      40      40       40      40      30      30      30      30
T      60  60      60      60      60      90      90      90      90       90      120     120     120     120     120     80      80      80       80      80      70      70      70      70
Z      20  20      20      20      20      30      30      30      30       30      20      20      20      20      20      10      10      10       10      10      12      12      12      12
  ;
parameter
co_tr(ti)/1 0.01,2 0.01,3 0.01,4 0.01,5 0.01,6 0.06,7 0.06,8 0.06,9 0.06,10 0.06,11 0.07,12 0.07,13 0.07,14 0.07,15 0.07,16 0.09,17 0.09,18 0.09,19 0.03,20 0.03,21 0.03,22 0.03,23 0.03,24 0.03/  ;


parameters cp_,cq_,pg_,qg_,pdo_,qdo_,pnl_,qnl_,Tr_,d_,dd;
*dd(r)=round(0.13*d(r));
loop((ti),
*Tr_(z,r,ti)=Tr(r,z);
D_(r,ti)=ceil(d(r)*co_tr(ti));
);

loop((ti),
cp_(ti)=co_price('T',ti);
cq_(ti)=co_price('T',ti)*0.15;
);
loop((i,ti),
pg_(i,ti)=co_p('T',ti)*pg(i);
qg_(i,ti)=co_p('T',ti)*qg(i);
pdo_(i,ti)=co_p('T',ti)*pdo(i);
qdo_(i,ti)=co_p('T',ti)*qdo(i);
);
loop((nl,ti),
pnl_(nl,ti)=co_p('T',ti)*pnl(nl);
qnl_(nl,ti)=co_p('T',ti)*qnl(nl);

)
parameter ha(k)/1  1,2  5,3   7,4 11,5 13,6 17,7 19,8 23,9 25,10 29,11 31/;

loop(k,
x_h(k,i,j)=ha(k)*branch(i,j,'x');
bb_h(k,i,j)=ha(k)*branch(i,j,'b');
)
variables
F
positive variable fr
positive variable landa
binary variable x
integer variable y
variables
istr
istim
inlrl
inliml
ir
iim
inlr
inlim
ilr
ilim
ilmg
ilrms
pst
*qst
sst
*pd(i)
*qd(i)
p
q
*THDv
vr
vim
istmg
istrms
vmg
vrms_v

ilr_loss
ilim_loss
ilr_sqr
ilim_sqr
**************
pp

ilim_mosbat
ilim_manfi
dilim


ilr_mosbat
ilr_manfi
dilr
;

Positive Variables
ilim_mosbat
ilim_manfi
ilim_sqr

ilr_mosbat
ilr_manfi
ilr_sqr;


*Integer Variables
*X
;

********************************

sos2 variables   lam1,lam2;
positive variables ast,bst,zst,f1,f2,gama ;

positive variables al,bl,zl,f1l,f2l,gamal;

positive variables av,bv,zv,f1v,f2v,gamav ,athd,bthd,zthd,f1thd,f2thd,gama_thd
equations
         obj_fun
         co1
         co2
         co3
         co4
         co5
         co6
         co7
         co8
         co9
         co10
         co11
         co12
         co13
*co14
*co15
         co16
         co17
         co18
         co19
         co20
         co21



         co33

         co40
         co41

cco1
cco2
cco3
cco4
cco5
cco6
cco7
*cco8
cco9
cco10
cco11
cco12
cco13
cco14
cco15
cco16
cco17
cco18
cco19
cco20
cco21
*cco22
cco23
cco24
cco25
cco26
cco27
cco28
cco29
cco30
cco31
cco32
cco33
cco34
cco35
cco36
cco37
cco38
cco39
cco40
cco41
cco42

cco44
cco45
cco46
cco47
cco48
cco49
cco50
cco51
cco52
cco53
cco54
cco55
cco56
cco57

cco59
cco60
cco61
cco62
cco63
cco64
cco65
cco66

*cccc1
cccc2
cs1

Eq2
Eq3
Eq4
Eq5
Eq6
Eq7
Eq8
eq9
eq10
eq11
eq12
eq13
eq14
eq15
eq16
eq17
;

parameters PEVs,S_PEVs,tlink,S_tlink,Ex_PEVs,S_Ex,Ex_t,S_Ex_t,
xx,yy,obj_funn,obj_fun_lc,c_instal,c_operation,c_travel,c_travel_LC,nnc,flow_link;

parameters ploss,ploss_t,vmg_p,vrms,THDv,qst,c_p,c_q,pp_loss,s_ilr;

****************************************G+BJ**********************************************
b_h(k,i,j)$line(i,j) = -x_h(k,i,j)/(sqr(branch(i,j,'r'))+sqr(x_h(k,i,j)));
g_h(k,i,j)$line(i,j) = branch(i,j,'r')/(sqr(branch(i,j,'r'))+sqr(x_h(k,i,j)));
ratio(line)=1;
************************************Ybus-Formulation***************************************************

y_h(k,i,j,'real')$(not sameas(i,j))=sum(Line(i,j)$branch(i,j,'st'), -1/ratio(i,j) * (g_h(k,i,j)*cos(branch(i,j,'an')) - b_h(k,i,j)*sin(branch(i,j,'an'))))
+ sum(Line(j,i)$branch(j,i,'st'), -1/ratio(j,i)* (g_h(k,j,i)*cos(-branch(j,i,'an')) - b_h(k,j,i)*sin(-branch(j,i,'an'))));


y_h(k,i,j,'imag')$(not sameas(i,j))=sum(Line(i,j)$branch(i,j,'st'), -1/ratio(i,j) * (b_h(k,i,j)*cos(branch(i,j,'an')) + g_h(k,i,j)*sin(branch(i,j,'an'))))
+ sum(Line(j,i)$branch(j,i,'st'), -1/ratio(j,i)* (b_h(k,j,i)*cos(-branch(j,i,'an')) + g_h(k,j,i)*sin(-branch(j,i,'an'))));


y_h(k,i,i,'real')= sum(j$branch(i,j,'st'), g_h(k,i,j)/sqr(ratio(i,j))) + sum(j$branch(j,i,'st'), g_h(k,j,i));



y_h(k,i,i,'imag')= sum(j$branch(i,j,'st'), 1/sqr(ratio(i,j)) * (b_h(k,i,j)+bb_h(k,i,j)/2)) + sum(j$branch(j,i,'st'),(b_h(k,j,i)+bb_h(k,j,i)/2));

***************************************objective_function***********************************************
obj_fun   ..F=e=sum(z,c_fix(z)*x(z)+cs_p(z)*y(z))+sum((a,z,ti),c_wait(z,a)*landa(z,a,ti))
+beta*(sum((z,r,ti),D_(r,ti)*fr(z,r,ti)*TR(r,z)))+sum((i,ti),cs(i)*Sbase*sst(i,ti))/24
+sum((ti),cp_(ti)*sbase*p('1',ti))
+sum((ti),cq_(ti)*sbase*q('1',ti))
+sum((i,j,ti)$ line(i,j) ,cp_(ti)*ilr_loss(i,j,ti)+ilim_loss(i,j,ti))
;
*****************************************TN_constraints*********************************************
co1(z,ti)                         ..sum(r,D_(r,ti)*fr(z,r,ti))=e=sum(a,landa(z,a,ti));
co2(z,a,ti)                      ..landa(z,a,ti)=l=deltan(z,a)*y(z);
co3(r,ti)                        ..sum(z,fr(z,r,ti))=e=1;
co4                           ..sum(z,x(z))=e=2 ;
co5(z)                        ..y(z)=l=10*x(z);
co6(z)                        ..y(z)=g=x(z);
co7(z,i,ti)$ node_bus(z,i)       ..sum(r,D_(r,ti)*fr(z,r,ti))=l=(psp(i)*deltat*y(z))/(cu*rng);

*****************************************************************************
co8(k,i ,ti)                   ..ir(k,i ,ti)=e=sum(j,y_h(k,i,j,'real')*vr(k,j ,ti)-y_h(k,i,j,'imag')*vim(k,j ,ti));
co9(k,i ,ti)                   ..iim(k,i ,ti)=e=sum(j,y_h(k,i,j,'real')*vim(k,j ,ti)+y_h(k,i,j,'imag')*vr(k,j ,ti));

co10(i ,ti)$lload(i)  ..ir('1',i ,ti)=e=(pg_(i ,ti)-(pdo_(i ,ti)+sum(nl$ nl_bus(i,nl),pnl_(nl ,ti))))*(2-vr('1',i ,ti))+(qg_(i ,ti)-(qdo_(i ,ti)+sum(nl$ nl_bus(i,nl),qnl_(nl ,ti))))*(vim('1',i ,ti))+istr('1',i ,ti);
co11(i ,ti)$lload(i)  ..iim('1',i ,ti)=e=(pg_(i ,ti)-(pdo_(i ,ti)+sum(nl$ nl_bus(i,nl),pnl_(nl ,ti))))*(vim('1',i ,ti))+(qg_(i ,ti)-(qdo_(i ,ti)+sum(nl$ nl_bus(i,nl),qnl_(nl ,ti))))*(-2+vr('1',i ,ti))+istim('1',i ,ti);

co12(i ,ti)$slack(i)                ..p(i ,ti)=e=ir('1',i ,ti);
co13(i ,ti)$slack(i)                ..q(i ,ti)=e=-iim('1',i ,ti);

co16(k,i ,ti)$(ord(k)<>1)           ..ir(k,i ,ti)=e=-istr(k,i ,ti)$lload(i)-inlr(k,i ,ti) ;
co17(k,i ,ti)$(ord(k)<>1)           ..iim(k,i ,ti)=e=-istim(k,i ,ti)$lload(i)-inlim(k,i ,ti);

co18(i,nl ,ti)$ (nl_bus(i,nl))      ..inlrl(i,nl ,ti)=e=pnl_(nl ,ti)*(2-vr('1',i ,ti))+qnl_(nl ,ti)*vim('1',i ,ti);
co19(i,nl ,ti)$(nl_bus(i,nl))       ..inliml(i,nl ,ti)=e=qnl_(nl ,ti)*(-2+vr('1',i ,ti))+pnl_(nl ,ti)*vim('1',i ,ti);

co20(k,i ,ti)$(ord(k)<>1)           ..inlr(k,i ,ti)=e=sum(nl,chr(nl,k)*inlrl(i,nl ,ti)$nl_bus(i,nl)-(chim(nl,k)*inliml(i,nl ,ti)$nl_bus(i,nl)));
co21(k,i ,ti)$(ord(k)<>1)           ..inlim(k,i ,ti)=e=sum(nl,chr(nl,k)*inliml(i,nl ,ti)$nl_bus(i,nl)+(chim(nl,k)*inlrl(i,nl ,ti)$nl_bus(i,nl)));

************************************************************************************************************************
cco37(k,i,d1 ,ti)               ..vr(k,i ,ti)*cos(360*k1(d1)/n)+vim(k,i ,ti)*sin(360*k1(d1)/n)=l=av(k,i ,ti);
cco38(i,d4 ,ti)                  ..av('1',i ,ti)*cos(360*k4(d4)/n)+av('2',i ,ti)*sin(360*k4(d4)/n)=l=bv('1',i ,ti);
cco39(i,d4 ,ti)                  ..av('3',i ,ti)*cos(360*k4(d4)/n)+av('4',i ,ti)*sin(360*k4(d4)/n)=l=bv('2',i ,ti);
cco40(i,d4 ,ti)                  ..av('5',i ,ti)*cos(360*k4(d4)/n)+av('6',i ,ti)*sin(360*k4(d4)/n)=l=bv('3',i ,ti);
cco41(i,d4 ,ti)                  ..av('7',i ,ti)*cos(360*k4(d4)/n)+av('8',i ,ti)*sin(360*k4(d4)/n)=l=bv('4',i ,ti);
cco42(i,d4 ,ti)                  ..av('9',i ,ti)*cos(360*k4(d4)/n)+av('10',i ,ti)*sin(360*k4(d4)/n)=l=bv('5',i ,ti);

cco44(i,d4 ,ti)                  ..av('11',i ,ti)*cos(360*k4(d4)/n)+f1v(i ,ti)*sin(360*k4(d4)/n)=l=1.05;
cco45(i,d4 ,ti)                  ..bv('1',i ,ti)*cos(360*k4(d4)/n)+bv('2',i ,ti)*sin(360*k4(d4)/n)=l=zv('1',i ,ti);
cco46(i,d4 ,ti)                  ..bv('3',i ,ti)*cos(360*k4(d4)/n)+bv('4',i ,ti)*sin(360*k4(d4)/n)=l=zv('2',i ,ti);
cco47(i,d4 ,ti)                  ..bv('5',i ,ti)*cos(360*k4(d4)/n)+f1v(i ,ti)*sin(360*k4(d4)/n)=l=zv('3',i ,ti);

cco48(i,d4 ,ti)                  ..zv('1',i ,ti)*cos(360*k4(d4)/n)+zv('2',i ,ti)*sin(360*k4(d4)/n)=l=gamav('1',i ,ti);

cco49(i,d4 ,ti)                  ..zv('3',i ,ti)*cos(360*k4(d4)/n)+f2v(i ,ti)*sin(360*k4(d4)/n)=l=gamav('2',i ,ti);

cco50(i,d4 ,ti)                  ..gamav('1',i ,ti)*cos(360*k4(d4)/n)+gamav('2',i ,ti)*sin(360*k4(d4)/n)=l=1.05;

**************************************************
cccc2(i ,ti)..                        (sum((ii,jj),lam2(i,ii,jj ,ti)*sqr(Voltage(ii,'Image')))+sum((ii,jj),lam1(i,ii,jj ,ti)*sqr(Voltage(ii,'Real'))))=g=0.903;
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cco51(i ,ti)                     ..vmg('1',i ,ti)=e=0.5*(sum((ii,jj),lam2(i,ii,jj ,ti)*sqr(Voltage(ii,'Image')))+sum((ii,jj),lam1(i,ii,jj ,ti)*sqr(Voltage(ii,'Real'))))+0.5;
*__________________________________________

cco52(k,i,d1 ,ti)$(ord(k)<> 1)    ..vr(k,i ,ti)*cos(360*k1(d1)/n)+vim(k,i ,ti)*sin(360*k1(d1)/n)=l=athd(k,i ,ti);

cco53(i,d4 ,ti)                  ..athd('1',i ,ti)*cos(360*k4(d4)/n)+athd('2',i ,ti)*sin(360*k4(d4)/n)=l=bthd('1',i ,ti);
cco54(i,d4 ,ti)                  ..athd('3',i ,ti)*cos(360*k4(d4)/n)+athd('4',i ,ti)*sin(360*k4(d4)/n)=l=bthd('2',i ,ti);
cco55(i,d4 ,ti)                  ..athd('5',i ,ti)*cos(360*k4(d4)/n)+athd('6',i ,ti)*sin(360*k4(d4)/n)=l=bthd('3',i ,ti);
cco56(i,d4 ,ti)                  ..athd('7',i ,ti)*cos(360*k4(d4)/n)+athd('8',i ,ti)*sin(360*k4(d4)/n)=l=bthd('4',i ,ti);
cco57(i,d4 ,ti)                  ..athd('9',i ,ti)*cos(360*k4(d4)/n)+athd('10',i ,ti)*sin(360*k4(d4)/n)=l=bthd('5',i ,ti);

cco59(i,d4 ,ti)                  ..athd('11',i ,ti)*cos(360*k4(d4)/n)+f1thd(i ,ti)*sin(360*k4(d4)/n)=l=0.05*vmg('1',i ,ti);
cco60(i,d4 ,ti)                  ..bthd('1',i ,ti)*cos(360*k4(d4)/n)+bthd('2',i ,ti)*sin(360*k4(d4)/n)=l=zthd('1',i ,ti);
cco61(i,d4 ,ti)                  ..bthd('3',i ,ti)*cos(360*k4(d4)/n)+bthd('4',i ,ti)*sin(360*k4(d4)/n)=l=zthd('2',i ,ti);
cco62(i,d4 ,ti)                  ..bthd('5',i ,ti)*cos(360*k4(d4)/n)+f1thd(i ,ti)*sin(360*k4(d4)/n)=l=zthd('3',i ,ti);
cco63(i,d4 ,ti)                  ..zthd('1',i ,ti)*cos(360*k4(d4)/n)+zthd('2',i ,ti)*sin(360*k4(d4)/n)=l=gama_thd('1',i ,ti);
cco64(i,d4 ,ti)                  ..zthd('3',i ,ti)*cos(360*k4(d4)/n)+f2thd(i ,ti)*sin(360*k4(d4)/n)=l=gama_thd('2',i ,ti);
cco65(i,d4 ,ti)                  ..gama_thd('1',i ,ti)*cos(360*k4(d4)/n)+gama_thd('2',i ,ti)*sin(360*k4(d4)/n)=l=0.05*vmg('1',i ,ti);
*________________________________________________________________
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cco66(k,i,d1 ,ti)$(ord(k)<> 1)    ..vr(k,i ,ti)*cos(360*k1(d1)/n)+vim(k,i ,ti)*sin(360*k1(d1)/n)=l=0.03*vmg('1',i ,ti);

*PEVFCS_constraints####

co33(z,i ,ti)$ (lload(i)and node_bus(z,i))        ..pst(i ,ti)=e=((sum(a,landa(z,a,ti))*cu*(0.001)*rng)/(deltat*bzde*sbase));

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cco30(i ,ti)$lload(i)                   ..pst(i ,ti)=e=sum((ii,jj),lam1(i,ii,jj ,ti)*Voltage(ii,'Real')*Current(jj,'Real'))+sum((ii,jj),lam2(i,ii,jj ,ti)*Voltage(ii,'Image')*Current(jj,'Image'));
cco31(i ,ti)                             ..vr('1',i ,ti)=e=sum((ii,jj),lam1(i,ii,jj ,ti)*Voltage(ii,'Real'));
cco32(i ,ti)                             ..istr('1',i ,ti)=e=sum((ii,jj),lam1(i,ii,jj ,ti)*Current(jj,'Real'));
cco33(i ,ti)                             ..sum((ii,jj),lam1(i,ii,jj ,ti))=e=1;
cco34(i ,ti)                             ..vim('1',i ,ti)=e=sum((ii,jj),lam2(i,ii,jj ,ti)*Voltage(ii,'Image'));
cco35(i ,ti)                             ..istim('1',i ,ti)=e=sum((ii,jj),lam2(i,ii,jj ,ti)*Current(jj,'Image'));
cco36(i ,ti)                             ..sum((ii,jj),lam2(i,ii,jj ,ti))=e=1;
*%__________________________________________________________________________
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cco1(z,i ,ti)$ (lload(i)and node_bus(z,i))                           ..sst(i ,ti)=l=M*x(z);

cco2(k,i,d1 ,ti)               ..istr(k,i ,ti)*cos(360*k1(d1)/n)+istim(k,i ,ti)*sin(360*k1(d1)/n)=l=ast(k,i ,ti);

cco3(i,d4 ,ti)                  ..ast('1',i ,ti)*cos(360*k4(d4)/n)+ast('2',i ,ti)*sin(360*k4(d4)/n)=l=bst('1',i ,ti);
cco4(i,d4 ,ti)                  ..ast('3',i ,ti)*cos(360*k4(d4)/n)+ast('4',i ,ti)*sin(360*k4(d4)/n)=l=bst('2',i ,ti);
cco5(i,d4 ,ti)                  ..ast('5',i ,ti)*cos(360*k4(d4)/n)+ast('6',i ,ti)*sin(360*k4(d4)/n)=l=bst('3',i ,ti);
cco6(i,d4 ,ti)                  ..ast('7',i ,ti)*cos(360*k4(d4)/n)+ast('8',i ,ti)*sin(360*k4(d4)/n)=l=bst('4',i ,ti);
cco7(i,d4 ,ti)                  ..ast('9',i ,ti)*cos(360*k4(d4)/n)+ast('10',i ,ti)*sin(360*k4(d4)/n)=l=bst('5',i ,ti);

cco9(i,d4 ,ti)                  ..ast('11',i ,ti)*cos(360*k4(d4)/n)+f1(i ,ti)*sin(360*k4(d4)/n)=l=sst(i ,ti)/1.05;

cco10(i,d4 ,ti)                  ..bst('1',i ,ti)*cos(360*k4(d4)/n)+bst('2',i ,ti)*sin(360*k4(d4)/n)=l=zst('1',i ,ti);
cco11(i,d4 ,ti)                  ..bst('3',i ,ti)*cos(360*k4(d4)/n)+bst('4',i ,ti)*sin(360*k4(d4)/n)=l=zst('2',i ,ti);
cco12(i,d4 ,ti)                  ..bst('5',i ,ti)*cos(360*k4(d4)/n)+f1(i ,ti)*sin(360*k4(d4)/n)=l=zst('3',i ,ti);
cco13(i,d4 ,ti)                  ..zst('1',i ,ti)*cos(360*k4(d4)/n)+zst('2',i ,ti)*sin(360*k4(d4)/n)=l=gama('1',i ,ti);
cco14(i,d4 ,ti)                  ..zst('3',i ,ti)*cos(360*k4(d4)/n)+f2(i ,ti)*sin(360*k4(d4)/n)=l=gama('2',i ,ti);
cco15(i,d4 ,ti)                  ..gama('1',i ,ti)*cos(360*k4(d4)/n)+gama('2',i ,ti)*sin(360*k4(d4)/n)=l=sst(i ,ti)/1.05;
*______________________________________________________________________________

co40(k,i,j ,ti)$ line(i,j)   ..ilr(k,i,j ,ti)=e=y_h(k,i,j,'real')*(vr(k,i ,ti)-vr(k,j ,ti))-y_h(k,i,j,'imag')*(vim(k,i ,ti)-vim(k,j ,ti));
co41(k,i,j ,ti)$ line(i,j)   ..ilim(k,i,j ,ti)=e=y_h(k,i,j,'real')*(vim(k,i ,ti)-vim(k,j ,ti))+y_h(k,i,j,'imag')*(vr(k,i ,ti)-vr(k,j ,ti));
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cco16(k,i,j,d1 ,ti)$ line(i,j)              ..ilr(k,i,j ,ti)*cos(360*k1(d1)/n)+ilim(k,i,j ,ti)*sin(360*k1(d1)/n)=l=al(k,i,j ,ti);

cco17(i,j,d4 ,ti)$ line(i,j)                  ..al('1',i,j ,ti)*cos(360*k4(d4)/n)+al('2',i,j ,ti)*sin(360*k4(d4)/n)=l=bl('1',i,j ,ti);
cco18(i,j,d4 ,ti)$ line(i,j)                 ..al('3',i,j ,ti)*cos(360*k4(d4)/n)+al('4',i,j ,ti)*sin(360*k4(d4)/n)=l=bl('2',i,j ,ti);
cco19(i,j,d4 ,ti)$ line(i,j)                  ..al('5',i,j ,ti)*cos(360*k4(d4)/n)+al('6',i,j ,ti)*sin(360*k4(d4)/n)=l=bl('3',i,j ,ti);
cco20(i,j,d4 ,ti)$ line(i,j)                  ..al('7',i,j ,ti)*cos(360*k4(d4)/n)+al('8',i,j ,ti)*sin(360*k4(d4)/n)=l=bl('4',i,j ,ti);
cco21(i,j,d4 ,ti)$ line(i,j)                  ..al('9',i,j ,ti)*cos(360*k4(d4)/n)+al('10',i,j ,ti)*sin(360*k4(d4)/n)=l=bl('5',i,j ,ti);

cco23(i,j,d4 ,ti)$ line(i,j)                  ..al('11',i,j ,ti)*cos(360*k4(d4)/n)+f1l(i,j ,ti)*sin(360*k4(d4)/n)=l=1.05;

cco24(i,j,d4 ,ti)$ line(i,j)                  ..bl('1',i,j ,ti)*cos(360*k4(d4)/n)+bl('2',i,j ,ti)*sin(360*k4(d4)/n)=l=zl('1',i,j ,ti);
cco25(i,j,d4 ,ti)$ line(i,j)                 ..bl('3',i,j ,ti)*cos(360*k4(d4)/n)+bl('4',i,j ,ti)*sin(360*k4(d4)/n)=l=zl('2',i,j ,ti);
cco26(i,j,d4 ,ti)$ line(i,j)                  ..bl('5',i,j ,ti)*cos(360*k4(d4)/n)+f1l(i,j ,ti)*sin(360*k4(d4)/n)=l=zl('3',i,j ,ti);
cco27(i,j,d4 ,ti)$ line(i,j)                  ..zl('1',i,j ,ti)*cos(360*k4(d4)/n)+zl('2',i,j ,ti)*sin(360*k4(d4)/n)=l=gamal('1',i,j ,ti);
cco28(i,j,d4 ,ti)$ line(i,j)                  ..zl('3',i,j ,ti)*cos(360*k4(d4)/n)+f2l(i,j ,ti)*sin(360*k4(d4)/n)=l=gamal('2',i,j ,ti);
cco29(i,j,d4 ,ti)$ line(i,j)                  ..gamal('1',i,j ,ti)*cos(360*k4(d4)/n)+gamal('2',i,j ,ti)*sin(360*k4(d4)/n)=l=1.05;

cs1(i ,ti)..            pst(i ,ti)=l=sst(i ,ti);
*****************************************************plosss_ilr******************************************************
Eq2(i,j,ti)$ line(i,j).. ilr_sqr('1',i,j ,ti)=e=sum(e,(2*ord(e)-1)*(ilr_max/ee)*dilr(e,i,j ,ti));
Eq3(i,j,ti)$ line(i,j).. ilr_mosbat(i,j,ti)-ilr_manfi(i,j,ti)=e=ilr('1',i,j ,ti);
Eq4(i,j,ti)$ line(i,j).. ilr_mosbat(i,j,ti)+ilr_manfi(i,j,ti)=e=sum(e,dilr(e,i,j,ti));
Eq5(e,i,j,ti)$ line(i,j).. dilr(e,i,j,ti)=g=0;
Eq6(e,i,j,ti)$ line(i,j).. dilr(e,i,j,ti)=l=dilr_max;
Eq7(i,j,ti)$ line(i,j).. ilr('1',i,j ,ti)=g=0;
Eq8(i,j,ti)$ line(i,j).. ilr('1',i,j ,ti)=l=ilr_max;

Eq9(i,j,ti)$ line(i,j).. ilr('1',i,j ,ti)*branch(i,j,'r')=e=ilr_loss(i,j,ti);

*****************************************************plosss_ilim******************************************************
Eq10(i,j,ti)$ line(i,j).. ilim_sqr('1',i,j ,ti)=e=sum(e,(2*ord(e)-1)*(ilim_max/ee)*dilim(e,i,j ,ti));
Eq11(i,j,ti)$ line(i,j).. ilim_mosbat(i,j,ti)-ilim_manfi(i,j,ti)=e=ilr('1',i,j ,ti);
Eq12(i,j,ti)$ line(i,j).. ilim_mosbat(i,j,ti)+ilim_manfi(i,j,ti)=e=sum(e,dilim(e,i,j,ti));
Eq13(e,i,j,ti)$ line(i,j).. dilim(e,i,j,ti)=g=0;
Eq14(e,i,j,ti)$ line(i,j).. dilim(e,i,j,ti)=l=dilim_max;
Eq15(i,j,ti)$ line(i,j).. ilr('1',i,j ,ti)=g=0;
Eq16(i,j,ti)$ line(i,j).. ilr('1',i,j ,ti)=l=ilim_max;

Eq17(i,j,ti)$ line(i,j).. ilim('1',i,j ,ti)*branch(i,j,'r')=e=ilim_loss(i,j,ti);


**********************************************************************************************************************************


        model dntnrecent /all/;
*option MINLP =lindo;
option MIP =cplex;
option optcr=0.1;
option optca=0.1;
         option Reslim=1000000000;
         solve dntnrecent using  MIP minimizing F;
**************DN*************
ploss(ti)=sum((k,i,j)$line(i,j),branch(i,j,'r')*( sqr(ilr.l(k,i,j,ti))+sqr(ilim.l(k,i,j,ti ))));
loop((i,k),
vmg_p(k,i,ti)=sqrt(sqr(vr.l(k,i,ti))+sqr(vim.l(k,i,ti))) ;);
loop( i,
vrms(i,ti)= sqrt(sum(k,sqr(vmg_p(k,i,ti)))););
loop(i,
THDv(i,ti)=(sqrt(sum(k$(ord(k)<> 1),sqr(vmg_p(k,i,ti)))))/(vmg.l('1',i,ti)););
loop((i,ti),
qst(i,ti)=(-vr.l('1',i,ti)*istim.l('1',i,ti))+(vim.l('1',i,ti)*istr.l('1',i,ti)); );
*******nonconsidering*********
    nnc(z,r,an)$ z_an(z,r,an)=1;
flow_link(z,r,an,ti)$ z_an(z,r,an)=x.l(z)*D_(r,ti)*fr.l(z,r,ti)*nnc(z,r,an);
PEVs(an,ti)=sum((z,r),flow_link(z,r,an,ti)$ z_an(z,r,an)) ;
Tlink(an,ti)=FFT(an)*(1+0.16*sqr(sqr(PEVs(an,ti)/(0.7*ca(an)))));
*************************************************
$ontext
*****considering*************************************
PEVs(an,ti)=vv.l(an,ti);
tlink(an,ti)=ttlink.l(an,ti);
$offtext
*********************************************
S_PEVs(ti)=sum(an,PEVs(an,ti));
S_tlink(ti)=sum(an,tlink(an,ti));
***************************************************
Ex_PEVs(an,ti)=max(0,PEVs(an,ti)-(round(0.7*ca(an))));
S_Ex(ti)=sum(an,Ex_PEVs(an,ti));

Ex_t(an,ti)=max(0,(tlink(an,ti)-fft(an)));
S_Ex_t(ti)=sum(an,Ex_t(an,ti));
*******************************************************

xx=sum(z,x.l(z));
yy=sum(z,y.l(z));
c_instal= sum(z,c_fix(z)*x.l(z)+cs_p(z)*y.l(z))+sum((i,ti),cs(i)*Sbase*sst.l(i,ti))/24   ;
c_operation=sum((a,z,ti),c_wait(z,a)*landa.l(z,a,ti));
c_travel=beta*(sum((z,r,an,ti)$ z_an(z,r,an),D_(r,ti)*fr.l(z,r,ti)*tlink(an,ti)));
c_travel_LC=beta*(sum((z,r,ti),D_(r,ti)*fr.l(z,r,ti)*Tr(r,z)));
obj_funn= c_instal+c_operation+c_travel;
obj_fun_lc=c_instal+c_operation+c_travel_LC;
c_p=-sum((ti),cp_(ti)*sbase*p.l('1',ti));
c_q=sum((ti),cq_(ti)*sbase*q.l('1',ti));
ploss_t=sum(ti,ploss(ti));
***************************************************
pp_loss=sum((i,j,ti)$ line(i,j) ,ilr_loss.l(i,j,ti)+ilim_loss.l(i,j,ti));
s_ilr(i,j,ti)=sqr(ilr.l('1',i,j ,ti));
*****************
display f.l,y.l,x.l,landa.l,fr.l,PEVs,S_PEVs,tlink,S_tlink,Ex_PEVs,S_Ex,Ex_t,S_Ex_t,
xx,yy,f.l,obj_funn,obj_fun_lc,c_instal,c_operation,c_travel,c_travel_LC,c_p,c_q,ploss_t
,p.l,q.l,ir.l,iim.l,inlr.l,inlim.l,istr.l,istim.l,sst.l,pst.l,qst,vr.l,vim.l,vrms,ploss,ploss_t,pp_loss,ilr_sqr.l,s_ilr,ilr.l;

execute_unload "results_56_no_c__with_ploss.gdx" VRMS
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=VRMS rng=classic1!'
execute_unload "results_56_no_c__with_ploss.gdx" ir
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx Var=ir rng=classic2!'
execute_unload "results_56_no_c__with_ploss.gdx" iim
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx Var=iim rng=classic3!'
execute_unload "results_56_no_c__with_ploss.gdx" THDv
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=THDv rng=classic4!'
execute_unload "results_56_no_c__with_ploss.gdx" pst
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx var=pst rng=classic5!'
execute_unload "results_56_no_c__with_ploss.gdx" qst
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=qst rng=classic6!'
execute_unload "results_56_no_c__with_ploss.gdx" sst
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx var=sst rng=classic7!'
execute_unload "results_56_no_c__with_ploss.gdx" ploss
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=ploss rng=classic8!'
**************************************************88
execute_unload "results_56_no_c__with_ploss.gdx" S_Ex
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=S_Ex rng=classic9!'
execute_unload "results_56_no_c__with_ploss.gdx" S_Ex_t
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=S_Ex_t rng=classic10!'

execute_unload "results_56_no_c__with_ploss.gdx" Ex_PEVs
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=Ex_PEVs rng=classic11!'
execute_unload "results_56_no_c__with_ploss.gdx" Ex_t
execute 'gdxxrw.exe results_56_no_c__with_ploss.gdx par=Ex_t rng=classic12!'

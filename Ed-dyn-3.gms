$TITLE This is a ramsey model with multiple sectors


$ontext

This model extends \DYN3.GMS to 
multiple sectors. As with DYN3.GMS 
the value returned to capital and labor 
in production is adjusted for consistency
with the observed investment level. A 
weighted least-sqrs. program is used to 
distribute the adjustment across the 
sectors.


Original Base Year Data (rectangular SAM)


       Y    X1  X2    I |     RA
----------------------------------------------------
PY   200            -20 |   -180
PX1 -120    120         |   
PX2  -80        80      |   
RENT            -20    -30      |     50
WAGE           -100    -50      |    150
Savings              20 |    -20

$offtext

Parameter
    Y0  Benchmark Gross Output,
    KS0 Benchmark Capital Supply,
    LS0 Benchmark Labor Supply,
    I0  Benchmark Savings,
    C0  Benchmark Consumption;

Y0  = 200;
KS0 =  50;
LS0 = 150;
I0  =  20;
C0  = 180;

Set J   Index on Production Sectors /1,2/;

Parameter
    X0(J)    Output by sector,
    LD(J)    Labor Demand by Sector,
    KD(J)    Capital Demand by Sector;

X0("1") =   120;
X0("2") =    80;
LD("1") =   100;
LD("2") =    50;
KD("1") =    20;
KD("2") =    30;

*----------Setup the Dynamics----------*
Set T       Time Periods    /2000*2100/
    TFIRST(T)   First Period    /2000/
    TLAST(T)    Terminal Period /2100/;

Scalars
    R   Rate of Interest    /0.05/
    G   Growth Rate     /0.02/,
    D   Depreciation Rate   /0.02/;

Parameter
    K0  Benchmark Capital Stock,
    RK0 Benchmark Return on a Unit of K
    PK0 Benchmark Price of a Unit of K
    QREF(T) Growth Path for Quantities
    PREF(T) Present Value Price Paths
    KSCAL0   Calibrated level of Value to Capital Services;

* 1) Given the interest rate Calculate PK0 and RK0

PK0 = 1/(1-r);
RK0 = (R-R*D+D)/(1-R);

* 2) Solve for the Initial Capital Stock

K0  = I0/(G+D);

* 3) Find the Calibrated Flow to Capital 

KSCAL0  =K0*RK0;

* 4) Adjust Social Accounts to match the new level
*     of capital services.  Here we employ a 
*     weighted-least-squares routine.  You may
*     want to use an alternative objective funct.?

POSITIVE VARIABLES 
    VK(J)   calibrated value of capital earnings
    VL(J)   calibrated value of labor earnings;

VARIABLE    OBJ objective function;

EQUATIONS
    VABAL(J)    value-added consistency by sector
    VKBAL       capital income balance
    OBJDEF      define the objective function;

PARAMETER   VA0(J)      Base year value added;

*   sectoral value-added

VA0(J) = LD(J) + KD(J);

*  Value of capital earnings consistent with the steady-state
*  growth path:

VABAL(J)..  VK(J) + VL(J)  =E= VA0(J);

VKBAL..     SUM(J, VK(J)) =E= KSCAL0;

OBJDEF..    OBJ =E= SUM((J), (1/KD(J))*SQR(VK(J) - KD(J)));

MODEL KBAL / VABAL, VKBAL, OBJDEF /;

KBAL.ITERLIM = 1000;
SOLVE KBAL USING NLP MINIMIZING OBJ;

*  Show what we did in the listing file.
DISPLAY "++++++++++++++++Raw factor Demands++++++++++++++++",LD, KD;

*   New Calibrated labor and capital demands:
LD(J) = VL.L(J);  
KD(J) = VK.L(J);

DISPLAY "+++++++++++++Calibrated Factor Demands++++++++++++",LD, KD;


* 5) Set the steady-state Reference Paths

QREF(T) =   (1+G)**(ORD(T)-1);
PREF(T) =   (1-r)**(ORD(T)-1);

Display K0,RK0,R,PK0,QREF,PREF;

Parameter TAX(T) Tax Rate on Capital Earnings;
TAX(T)  =  0;

$ONTEXT

$MODEL:DYN4

$Sectors:
    Y(T)    ! Macro Output (transitory utility)
    X(J,T)  ! Production
    I(T)    ! Investment
    K(T)    ! Capital Stock
    C(T)    ! Consumption Index

$Commodities:
    PY(T)   ! Price index on macro output
    PX(J,T) ! Price index on sector output
    PC(T)   ! Price index on consumption
    RK(T)   ! Present Value Return to capital
    PL(T)   ! Present Value Wage
    PK(T)   ! Price index on Capital
    PKT ! Price of Terminal Capital

$Consumers:
    RA  ! Representative agent

$Auxiliary:
    TCAP    ! Terminal Capital Demand

$PROD:Y(T) s:1
    O:PY(T)   Q:Y0
    I:PX(J,T) Q:X0(J)

$PROD:X(J,T) s:1
    O:PX(J,T) Q:X0(J)
    I:RK(T)   Q:KD(J)
    I:PL(T)   Q:LD(J)

$PROD:I(T)
    O:PKT$TLAST(T)  Q:I0
    O:PK(T+1)   Q:I0
    I:PY(T)     Q:I0    A:RA T:TAX(T)

$PROD:K(T)
    O:PKT$TLAST(T)  Q:(K0*(1-D))
    O:PK(T+1)   Q:(K0*(1-D))
    O:RK(T)     Q:KSCAL0
    I:PK(T)     Q:K0

$PROD:C(T) 
    O:PC(T)     Q:C0
    I:PY(T)     Q:C0

$Demand:RA s:0.5
    D:PC(T)     Q:(C0*QREF(T))  P:PREF(T)
    E:PK(TFIRST)    Q:K0
    E:PL(T)     Q:(SUM(j,LD(j))*QREF(T))
    E:PKT       Q:(-1)      R:TCAP

$Constraint:TCAP
    SUM(T$TLAST(T+1),C(T)*I(T+1) - I(T)*C(T+1)) =E= 0;

$offtext
$sysinclude mpsgeset dyn4

* Set Steady State Level Values:
Y.L(T)  = QREF(T);
X.L(J,T)= QREF(T);
I.L(T)  = QREF(T);
K.L(T)  = QREF(T);
C.L(T)  = QREF(T);

PY.L(T) = PREF(T);
PX.L(J,T)= PREF(T);
PC.L(T) = PREF(T);
RK.L(T) = PREF(T);
PL.L(T) = PREF(T);
PK.L(T) = PREF(T)*PK0;


PKT.L   = (1-R)*SUM(TLAST, PK.L(TLAST));
TCAP.L  = SUM(TLAST, I0*QREF(TLAST)+K0*(1-D)*QREF(TLAST));

dyn4.ITERLIM = 0;
$include dyn4.gen
solve dyn4 using mcp;

* Run a simple experiment: 10% tax on investment from 2010 on.

SET 
    SHOCK(T) Time periods under the shock /2010*2100/;

TAX(SHOCK) = 0.1;

dyn4.ITERLIM = 100;
$include dyn4.gen
solve dyn4 using mcp;
$exit

Parameter MACRO(T,*)    "% Change in Macro Variables",
      OUTPUT(T,*)   "% Change in output";
MACRO(T,"INVEST")   = 100*(I.L(T)/QREF(T)-1);
MACRO(T,"CONS")     = 100*(C.L(T)/QREF(T)-1);
MACRO(T,"CAPITAL")  = 100*(K.L(T)/QREF(T)-1);
MACRO(T,"OUTPUT")   = 100*(Y.L(T)/QREF(T)-1);

OUTPUT(T,"Aggregate")   = 100*(Y.L(T)/QREF(T)-1);
OUTPUT(T,"Sector1") = 100*(X.L("1",T)/QREF(T)-1);
OUTPUT(T,"Sector2") = 100*(X.L("2",T)/QREF(T)-1);

$libinclude gnuplot macro
$libinclude gnuplot output
Display macro,output;

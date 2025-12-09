%let pgm=utl-altair-slc-using-unit-square-mapping-calculate-max-pKa-slope-of-titration-curve-for-diprotic-acid;

%stop_submission;

Altair slc using unit square mapping calculate max pKa slope of titration curve for diprotic acid

see high res graph
Comparison of Rogers Fit to Best Stackoverflow fit (plot created by altair slc sglot in d:/png/nlin_ploots.png)
https://github.com/rogerjdeangelis/utl-altair-slc-using-unit-square-calculate-pKa-maximum-slope-from-titration-curve-for-diprotic-acid/blob/main/nlin_plots.png

Too long to poste on listserv

see github
https://github.com/rogerjdeangelis/utl-altair-slc-using-unit-square-calculate-pKa-maximum-slope-from-titration-curve-for-diprotic-acid

A new general framework for fitting non linear equations with transcendental and beta cdf/pdf functions on the unit square
The beta pdfs and cdfs are by far the most flexible formulations for fitting non-linear  functions over a unit square.
Note the beta has a natural relation to the polynomials, the expansion of the beta, x^ (1-x)^  is a polynomial.

CONTENTS
    1 Roger Solution
       map points to unit square
       add weight around near vertical segment (regression is not the best soultion to fit vertical lines?)
       run slc non-linear regression (proc nlin) on unit square mapping
       plot solution in unit square
       map unit square back to orginal units
    2  Stackoverflow solution
    3  Comparison of Rogers solution to stackoverflow solution
    4  Calculate pKa equilibrium

stackoverflow
https://stackoverflow.com/questions/78915512/using-python-to-calculate-pka-from-titration-curve-for-diprotic-acid

 Regression cannot easily handle an almost perfectly vertical line segment.
 Regression may not be appropriate for this problem. However, this solution
 provides the methodology to handle other complex non-linear fits.

 What is remarkable about this fit is that it only needs 3 coefficients.
 I have seen models using a knot that do not seem to fit that much better with 8 parameters.
 A knot provides for a knot and two equations and two sets of parameters

A general framework for non-linear fits?

1 The unit square refers to transforming your data so that all
  values fall within a square with sides of length 1 (from 0 to 1 on both axes). This normalization
  technique can make certain mathematical operations and comparisons easier.
  The unit square has special significance ii complex numbers.

2 Transcendental functions whose Y aor X domains often map to the unit square and can adapt to many data patterns.
  These functions are often orthogonal and independent of each other, which can be useful in modeling complex
  relationships.

3 Transcendental functions include sine, cosine, and tangent functions. They're periodic and can model
  cyclical patterns in data. Which are independent, orthogonal and partition the mean square error
  exclusively. They can naturally map into a unit square,

4 The simple beta, cdf, and pdf, have only two parameters and is well defined in the unit square.
  The distribution function (CDF) and probability density function (PDF) can model various
  shapes. The most flexible of all cdfs and pdfs?

  The combination of these techniques
  and functions allow for modeling many different types of relationships and patterns in data.


RELATED UNIT QSQUARE SOLUTIONS

https://github.com/rogerjdeangelis/utl-general-framework-for-fitting-non-linear-equation-with-transcendental-functions-on-a-unit-square
https://github.com/rogerjdeangelis/utl-better-fit-with-two-paramater-beta-then-four-than-four-polynomial-on-unit-square
https://github.com/rogerjdeangelis/utl-piecewise-linear-fit-compared-to-unit-square-beta-cdf-non-linear-fit-two-vs-four-parameters
https://github.com/rogerjdeangelis/utl-using-beta-cdf-and-pdf-unitsquare-to-fit-non-linear-equations-and-simple-polynomials
https://github.com/rogerjdeangelis/utl-using-the-unitsquare-framework-and-cumulative-beta-to-fit-a-non-linear-equation

/******************************************************************************************************************/
/*                                                                                                                */
/*   POSTED SOLUTION ON STACKOVERFLOW COMPARED TO MY SOLUTION                                                     */
/*                                                                                                                */
/*                                                                                                                */
/*       --+---------+---------+---------+---------+---------+---------+---------+---------+---------+--          */
/*       |                                                                                             |          */
/*       |  OPS FIT ='o'                    ROGER = 'r'                   ORIGINAL DATA='.'            |          */
/*       |                                                                                             |          */
/*       |  a=8.8179917838                  A =  307.4                                                 |          */
/*       |  b=0.9573331691                  B =  500.4                                                 |          */
/*       |  c=341.0840722484                C = 0.0422                                                 |          */
/*       |  d=2.3689668227                                                                             |          */
/*       |  e=1 ;                                                                                      |          */
/*       |                                                                                             |          */
/*       |  py = A/(1 + B**(xs - C)) + D    Ymapped =sqrt(x)*CDF('BETA',x,a,b) + c*sin(1.57*x) [0,1]   |          */
/*       |                                                                                             |          */
/*       |                                  FINAL EQUATION IN ORIGINAL UNITS                           |          */
/*       |                                                                                             |          */
/*       |                                   XT=x *867;                                                |          */
/*       |                                   YT =2.22 + (sqrt(x)*CDF('BETA',x,307.3,500.4)             |          */
/*       |                                          + 0.04218*sin(1.57*x))*(12.49-2.22); [0,900]       |          */
/*       |                                                                                             |          */
/*       |                                                                                             |          */
/*       | 0        100       200       300       400  XT   500       600       700       800       900|          */
/*       --+---------+---------+---------+---------+---------+---------+---------+---------+---------+--          */
/*    YT |                                                                                 rrrr ....   |          */
/*       |                                                                               rr....        |          */
/*  12.2 +  'o'=ops                                                                 ........           + 12.2     */
/*       |  'r'=rogers fit                                                     .....rr                 |          */
/*       |  '.'=orginal data                                                ....rrr                    |          */
/*       |                                                                ...rrr                       |          */
/*       |                                                              ...rr                          |          */
/*  11.2 +                                              oooooooooooooo...oooooooooooooooooooooooooo    + 11.2     */
/*       |                                           oooo           r..     stackoverflow fit          |          */
/*       |                          'o'=ops fit-->  oo           rr..                                  |          */
/*       |                                         oo          rr...                                   |          */
/*       |                                        oo        rrr...                                     |          */
/*  10.2 +                                        o       rr...                                        + 10.2     */
/*       |                                       oo     r....   <- '.'  original data                  |          */
/*       |                                       o   rr...                                             |          */
/*       |                                       o r....                                               |          */
/*       |                                      or...                                                  |          */
/*   9.2 +                                      o..                                                    +  9.2     */
/*       |                                     ...                                                     |          */
/*       |                 'r'=rogers fit --> r.o                                                      |          */
/*       |                                    .o                                                       |          */
/*       |                                    .o                                                       |          */
/*   8.2 +                                   .ro                                                       +  8.2     */
/*       |                                   .oo                                                       |          */
/*       |                                   .o                                                        |          */
/*       |                                   ro                                                        |          */
/*       |                                   .o                                                        |          */
/*   7.2 +                                   .o                                                        +  7.2     */
/*       |                                   o                                                         |          */
/*       |                                  .o                                                         |          */
/*       |                                   o                                                         |          */
/*       |                                  ro                                                         |          */
/*   6.2 +                                  .o                                                         +  6.2     */
/*       |   pKa equlibrium point           o  XT=330 YT=5.73                                          |          */
/*       |----------------------------------+----------------------------------------------------------|          */
/*       |                                  |                                                          |          */
/*       |                                  o                                                          |          */
/*   5.2 +                                 o|                                                          +  5.2     */
/*       |                                 or                                                          |          */
/*       |                                 or                                                          |          */
/*       |                                 o.                                                          |          */
/*       |                                or.                                                          |          */
/*   4.2 +                                or.                                                          +  4.2     */
/*       |                                or.                                                          |          */
/*       |                               oo..  <- '.' Original data                                    |          */
/*       |                               o .|                                                          |          */
/*       |                              o..r|                                                          |          */
/*   3.2 +                              ..r |                                                          +  3.2     */
/*       |                           .... r |                                                          |          */
/*       |                          ..oo  r |                                                          |          */
/*       |                 .......  ooo  rr |                                                          |          */
/*       | ooooo............oooooo  rrrrrr  |                                                          |          */
/*   2.2 + ......rrrrr                      |                                                          +  2.2     */
/*       |                                  |                                                          |          */
/*       --+---------+---------+---------+---------+---------+---------+---------+---------+---------+--          */
/*         0        100       200       300       400       500       600       700       800       900           */
/*                                                     XT                                                         */
/******************************************************************************************************************/
/*                   _
(_)_ __  _ __  _   _| |_
| | `_ \| `_ \| | | | __|
| | | | | |_) | |_| | |_
|_|_| |_| .__/ \__,_|\__|
        |_|
*/

/*--- CREAR OUT ALTERNATE PERMANENT WORK LIBRARY ---*/

proc datasets lib=workx kill;   /*--- clear d:/wpswrkx ---*/
run;quit;

Download into "d:/wpswrkx" ie d:/wpswrkx/rwant.sas7bdat. This is a permanent alternate work library;
https://github.com/rogerjdeangelis/utl-altair-slc-using-unit-square-calculate-pKa-maximum-slope-from-titration-curve-for-diprotic-acid/blob/main/rwant.sas7bdat

d:/wpswrkx/rwant.sas7bdat
workx.RWANT total obs=868 09DEC2025:10:27:40

 ROWNAMES    XS      YS

     1        0    2.2200
     2        1    2.2215
     3        2    2.2230
     4        3    2.2245
     5        4    2.2260
 ...
   864      863    12.482
   865      864    12.484
   866      865    12.486
   867      866    12.488
   868      867    12.490

/*--- copy sd1

options ls=90 ps=44;
proc plot data= workx.rwant;
 plot ys*xs='*' /box;
run;quit;
options ps=65 ls=171;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                              XS                                                                        */
/*        0        110       220       330       440       550       660       770       880                              */
/*     ---+---------+---------+---------+---------+---------+---------+---------+---------+---                            */
/*     |                                                                                     |                            */
/*  14 +                                                                                     + 14                         */
/*     |   pKa curve                                                                         |                            */
/*     |                                                                                     |                            */
/*     |   SD1.RWANT total obs=868                                                           |                            */
/*     |                                                                       ***********   |                            */
/*  12 +     XS      YS                                               *********              + 12                         */
/*     |                                                          *****                      |                            */
/*     |      0    2.2200                                      ****                          |                            */
/*     |      1    2.2215                                   ****                             |                            */
/*     |      2    2.2230                               *****                                |                            */
/*  10 +     ..                                    ******                                    + 10                         */
/*     |     865   12.484                     ******                                         |                            */
/*     |     866   12.486                   ***                                              |                            */
/*  YS |     867   12.488                 ***                                                | YS                         */
/*     |     868   12.490                **                                                  |                            */
/*   8 +                                 *                                                   + 8                          */
/*     |                                 *                                                   |                            */
/*     |                                *                                                    |                            */
/*     |                                *  pKa of an acid from titration                     |                            */
/*     |                                *                                                    |                            */
/*   6 +                                *  Difficult to fit a very near vertical line        + 6                          */
/*     |                                   with non-linear or linear regression?             |                            */
/*     |                                *                                                    |                            */
/*     |                                   I don't see how a physical process can            |                            */
/*     |                                *  instantly have two y states?                      |                            */
/*   4 +                               **                                                    + 4                          */
/*     |                              **                                                     |                            */
/*     |                           **** |                                                    |                            */
/*     |                    ********    |  PKa equivalence point                             |                            */
/*     |  *******************           |  ~330                                              |                            */
/*   2 +                                |                                                    + 2                          */
/*     |                                |                                                    |                            */
/*     ---+---------+---------+---------+---------+---------+---------+---------+---------+---                            */
/*        0        110       220       330       440       550       660       770       880                              */
/*                                                                                                                        */
/**************************************************************************************************************************/
/*                                         _       _   _
/ |  _ __ ___   __ _  ___ _ __   ___  ___ | |_   _| |_(_) ___  _ __
| | | `__/ _ \ / _` |/ _ \ `__| / __|/ _ \| | | | | __| |/ _ \| `_ \
| | | | | (_) | (_| |  __/ |    \__ \ (_) | | |_| | |_| | (_) | | | |
|_| |_|  \___/ \__, |\___|_|    |___/\___/|_|\__,_|\__|_|\___/|_| |_|
                |___/    _                      _ _
 _ __ ___   __ _ _ __   | |_ ___    _   _ _ __ (_) |_   ___  __ _ _   _  __ _ _ __ ___
| `_ ` _ \ / _` | `_ \  | __/ _ \  | | | | `_ \| | __| / __|/ _` | | | |/ _` | `__/ _ \
| | | | | | (_| | |_) | | || (_) | | |_| | | | | | |_  \__ \ (_| | |_| | (_| | | |  __/
|_| |_| |_|\__,_| .__/   \__\___/   \__,_|_| |_|_|\__| |___/\__, |\__,_|\__,_|_|  \___|
                |_|                                            |_|
*/

proc sql;
  create
     table workx.nrmbfr as
  select
     xs as xorg
    ,ys as yorg
    ,(-min(xs)+xs)/(max(xs)-min(xs)) as xs
    ,(-min(ys)+ys)/(max(ys)-min(ys)) as ys
  from
     workx.rwant
;quit;


options ls=90 ps=44;
proc plot data= workx.nrmbfr;
 plot ys*xs='*' /box haxis=0 to 1 by .2;
run;quit;
options ps=65 ls=171;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*                             Plot of YS*XS.  Symbol used is '*'.                                                        */
/*                                                                                                                        */
/*      --+---------------+---------------+---------------+---------------+---------------+---                            */
/*   YS |                                                                                    |                            */
/*      |  MAPPED TO UNIT SQUARE                                                             |                            */
/*  1.0 +                                                                          ********  +                            */
/*      |  WORK.NRMBFR total obs=868                                      **********         |                            */
/*      |                                                            ******                  |                            */
/*      |  ROW       XS          YS                            ****                          |                            */
/*      |                        ****                                                        |                            */
/*      |    0    0.000000    .0000000                    ***                                |                            */
/*  0.8 +    1    0.001153    .0001461                *****                                  +                            */
/*      |    2    0.002307    .0002921            *****                                      |                            */
/*      |    3    0.003460    .0004382           *****                                       |                            */
/*      |                                     ****                                           |                            */
/*      |  865    0.99769      0.99961      ***                                              |                            */
/*      |  866    0.99885      0.99981    **                                                 |                            */
/*  0.6 +  867    1.00000      1.00000    *                                                  +                            */
/*      |                                **                                                  |                            */
/*      |                                *                                                   |                            */
/*      |                                *                                                   |                            */
/*      |                                                                                    |                            */
/*      |                                *                                                   |                            */
/*  0.4 +                                *                                                   +                            */
/*      |                                                                                    |                            */
/*      |                                *                                                   |                            */
/*      |                                *                                                   |                            */
/*      |                                                                                    |                            */
/*      |                               *                                                    |                            */
/*  0.2 +                               *                                                    +                            */
/*      |                               *                                                    |                            */
/*      |                              **                                                    |                            */
/*      |                           ****                                                     |                            */
/*      |                      ******                                                        |                            */
/*      |           ************                                                             |                            */
/*  0.0 + ***********                                                                        +                            */
/*      |                                                                                    |                            */
/*      --+---------------+---------------+---------------+---------------+---------------+---                            */
/*       0.0             0.2             0.4             0.6             0.8             1.0                              */
/*                                                                                                                        */
/**************************************************************************************************************************/


/*                                            __ _ _   _   _
 _ __  _ __ ___   ___ ___  ___ ___           / _(_) |_| |_(_)_ __   __ _
| `_ \| `__/ _ \ / __/ _ \/ __/ __|  _____  | |_| | __| __| | `_ \ / _` |
| |_) | | | (_) | (_|  __/\__ \__ \ |_____| |  _| | |_| |_| | | | | (_| |
| .__/|_|  \___/ \___\___||___/___/         |_| |_|\__|\__|_|_| |_|\__, |
|_|_     _                _       _     _                           |___/
  __ _  __| | __| | __      _____(_) __ _| |__ | |_ ___
 / _` |/ _` |/ _` | \ \ /\ / / _ \ |/ _` | `_ \| __/ __|
| (_| | (_| | (_| |  \ V  V /  __/ | (_| | | | | |_\__ \
 \__,_|\__,_|\__,_|   \_/\_/ \___|_|\__, |_| |_|\__|___/
                                    |___/
*/
                           /*--- add weight around the near vertical segment ---*/
data workx.addweights;
 set workx.nrmbfr;
  weights=1;
  if .35 <xs < .38 then weights=10;
  if .39 <xs < .45 then weights=10;
run;quit;


/*__ _ _
 / _(_) |_
| |_| | __|
|  _| | |_
|_| |_|\__|

*/

ods trace on;
ods output parameterestimates=workx.est;
proc nlin data=workx.addweights method=MARQUARDT CONVERGE=0.1;;
  parms A = 8
        B = 5
        c = 2
        ;
  model ys = sqrt(xs)*CDF('BETA', xs, a, b) + c*sin(1.57*xs) ;
  bounds A > 0, B > 0;
  output out=workx.results p=predicted r=residual;
run;quit;
ods trace off;

proc print data=workx.est;
format _numeric_ 8.5;
run;quit;

options ls=90 ps=44;
proc plot data=workx.results;
 plot residual*xs='*' /box;
 plot   ys*xs='*' predicted*xs='P' ys*xs='*' / overlay box
         haxis=0 to 1 by .1
         vaxis=0 to 1 by .1;
 plot predicted*ys='*' /box;
run;quit;

/*           _               _                                         _
  ___  _   _| |_ _ __  _   _| |_  _ __ ___   __ _ _ __  _ __   ___  __| |
 / _ \| | | | __| `_ \| | | | __|| `_ ` _ \ / _` | `_ \| `_ \ / _ \/ _` |
| (_) | |_| | |_| |_) | |_| | |_ | | | | | | (_| | |_) | |_) |  __/ (_| |
 \___/ \__,_|\__| .__/ \__,_|\__||_| |_| |_|\__,_| .__/| .__/ \___|\__,_|
                |_|                              |_|   |_|
*/

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  Only 4 parameters                                                                                                     */
/*                                                                                                                        */
/*  MODEL  YS = YS = SQRT(XS) * CDF('BETA', XS, A, B) + C*SIN(1.57*XS)                                                    */
/*                                                                                                                        */
/*                                    Sum of        Mean               Approx                                             */
/*  Source                    DF     Squares      Square    F Value    Pr > F                                             */
/*                                                                                                                        */
/*  Model                      3       395.9       132.0     210152    <.0001                                             */
/*  Error                    865      0.5431    0.000628                                                                  */
/*  Uncorrected Total        868       396.4                                                                              */
/*                                                                                                                        */
/*                                         Approximate 95%                                                                */
/*                             Approx         Confidence                      Approx                                      */
/* PARAMETER    ESTIMATE    Std Error           Limits           t Value    Pr > |t|                                      */
/*                                                                                                                        */
/* A               309.6      17.4391       275.4       343.8      17.75      <.0001                                      */
/* B               504.0      28.3750       448.3       559.7      17.76      <.0001                                      */
/* c              0.0422      0.00123      0.0398      0.0446      34.39      <.0001                                      */
/*                                                                                                                        */
/*  ys=sqrt  0.0422      0.00123(xs)*CDF('BETA',xs,a,b)+c*sin(1.57*xs)                                                    */
/*                                                                                                                        */
/*  PREDICTED AND OBSERVED DATA                                                                                           */
/*                                                                                                                        */
/*                                                                                                                        */
/*              0.0    0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1.0                                 */
/*            ---+------+------+------+------+------+------+------+------+------+------+---                               */
/*            |                                                                           |                               */
/*            |  PREDICTED AND OBSERVED DATA                                              |                               */
/*            |                                                                           |                               */
/*            |                                                                           |                               */
/*            |  y =sqrt(x)*CDF('BETA',x,a,b) + c*in(1.57*x)                              |                               */
/*            |                                                                           |                               */
/*            |  * = Original mapped data                                                 |                               */
/*            |  p = Original and mapped poins (most of the points)                       |                               */
/*  PREDICTED |                                                                     PPPP  |                               */
/*        1.0 +                                                                PPPPP****  + 1.0                           */
/*            |                                                          **PPPPP**        |                               */
/*            |                                                      **PPPPP              |                               */
/*        0.9 +                                                   *PPPPP                  + 0.9                           */
/*            |                                                PPPPP                      |                               */
/*            |                                             PPPP*                         |                               */
/*        0.8 +                                          PPPP*                            + 0.8                           */
/*            |                                      PPPP**                               |                               */
/*            |                                   PPPP**                                  |                               */
/*        0.7 +                                PPPP*                                      + 0.7                           */
/*            |                               PP*                                         |                               */
/*            |                              PP                                           |                               */
/*        0.6 +                              P                                            + 0.6                           */
/*            |                             *P                                            |                               */
/*            |                             *P                                            |                               */
/*        0.5 +                             PP                                            + 0.5                           */
/*            |                             P                                             |                               */
/*            |                             P                                             |                               */
/*        0.4 +                             P                                             + 0.4                           */
/*            |                             P                                             |                               */
/*            |                             P                                             |                               */
/*        0.3 +                            PP                                             + 0.3                           */
/*            |                            P*                                             |                               */
/*            |                            P*                                             |                               */
/*        0.2 +                            P*                                             + 0.2                           */
/*            |                            P                                              |                               */
/*            |                           PP                                              |                               */
/*        0.1 +                         **P                                               + 0.1                           */
/*            |                     ***** P                                               |                               */
/*            |           *********PPPPPPPP                                               |                               */
/*        0.0 +  PPPPPPPPPPPPPPPPPPP                                                      + 0.0                           */
/*            ---+------+------+------+------+------+------+------+------+------+------+---                               */
/*              0.0    0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1.0                                 */
/*                                                                                                                        */
/*                                                 XS                                                                     */
/*  RESUDUALS                                                                                                             */
/*                                                                                                                        */
/*                                                                                                                        */
/*                          Plot of RESIDUAL*XS.  Symbol used is '*'.                                                     */
/*                                                                                                                        */
/*                                                                                                                        */
/*            0.0000   0.1244   0.2488   0.3733   0.4977   0.6221   0.7465   0.8710   0.9954                              */
/*            ---+--------+--------+--------+--------+--------+--------+--------+--------+---                             */
/*   RESIDUAL |                                                                             |                             */
/*       0.10 +                                                                             +  0.10                       */
/*            |  RESUDUALS VS OBSERVED X                                                    |                             */
/*            |                                                                             |                             */
/*            |                              *                                              |                             */
/*            |                          **  *                                              |                             */
/*            |                         **** **                                             |                             */
/*       0.05 +                        **  *  *                                             +  0.05                       */
/*            |                      ***   *  *                                             |                             */
/*            |                    ***     *  *                                             |                             */
/*            |                  ***       * **                                             |                             */
/*            |              ****          *  *                     *********               |                             */
/*            |       *******              *  *                    **       ****            |                             */
/*       0.00 +  ******                    *   *                 ***           **           +  0.00                       */
/*            |                            *   * ******         **               ***        |                             */
/*            |                             *  ***    **********                   ***      |                             */
/*            |                             *                                        ***    |                             */
/*            |                             **                                         **   |                             */
/*            |                             *                                           **  |                             */
/*      -0.05 +                                                                             + -0.05                       */
/*            |                             **                                              |                             */
/*            |                             *                                               |                             */
/*            |                             *                                               |                             */
/*            |                              *  pKa equivalece point (issue)                |                             */
/*            |                             *                                               |                             */
/*      -0.10 +                             **                                              + -0.10                       */
/*            |                                                                             |                             */
/*            |                                                                             |                             */
/*            ---+--------+--------+--------+--------+--------+--------+--------+--------+---                             */
/*            0.0000   0.1244   0.2488   0.3733   0.4977   0.6221   0.7465   0.8710   0.9954                              */
/*                                                                                                                        */
/*                                                  XS                                                                    */
/*                                                                                                                        */
/*   PREDICTED VS OBSERVED DEPENDANT VARIABLE(Y)                                                                          */
/*                                                                                                                        */
/*             ---+------------+------------+------------+------------+------------+---                                   */
/*   PREDICTED |                                                                      |                                   */
/*         1.2 +  PREDICTED VS OBSERVED DEPENDANT VARIABLE(Y)                         +                                   */
/*             |                                                                      |                                   */
/*             |                                                                      |                                   */
/*             |                                                                      |                                   */
/*             |                                                                   *  |                                   */
/*         1.0 +                                                                 ***  +                                   */
/*             |                                                                **    |                                   */
/*             |                                                              ***     |                                   */
/*             |                                                          *****       |                                   */
/*             |                                                      *****           |                                   */
/*         0.8 +                                                    ***               +                                   */
/*             |                                                 ****                 |                                   */
/*             |                                               ***                    |                                   */
/*             |                                            ****                      |                                   */
/*             |                                           **                         |                                   */
/*         0.6 +                                          **                          +                                   */
/*             |                                        **                            |                                   */
/*             |                                       **                             |                                   */
/*             |                                      *                               |                                   */
/*             |                                   ***                                |                                   */
/*         0.4 +                             *    *                                   +                                   */
/*             |                   *  *  *                                            |                                   */
/*             |                **                                                    |                                   */
/*             |             ***                                                      |                                   */
/*             |            *                                                         |                                   */
/*         0.2 +           **                                                         +                                   */
/*             |           *                                                          |                                   */
/*             |          **                                                          |                                   */
/*             |         **                                                           |                                   */
/*             |       ***                                                            |                                   */
/*         0.0 +  ******                                                              +                                   */
/*             |                                                                      |                                   */
/*             ---+------------+------------+------------+------------+------------+---                                   */
/*               0.0          0.2          0.4          0.6          0.8          1.0                                     */
/*                                                                                                                        */
/*                                                YS                                                                      */
/*                                                                                                                        */
/**************************************************************************************************************************/
/*
 _                        __                       _                _
| |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___  | |__   __ _  ___| | __
| __| `__/ _` | `_ \/ __| |_ / _ \| `__| `_ ` _ \ | `_ \ / _` |/ __| |/ /
| |_| | | (_| | | | \__ \  _| (_) | |  | | | | | || |_) | (_| | (__|   <
 \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_||_.__/ \__,_|\___|_|\_\


data workx.back;
  set workx.results;
  xt=xs*867;
  yt = 2.22 + (12.49-2.22)*ys;
  pt = 2.22 + (12.49-2.22)*predicted ;
  keep xt yt pt ox oy;
run;quit;

options ls=90 ps=44;
proc plot data= workx.back;
 plot yt*xt='*' /box vaxis=2.20 to 12.5 by 1;
run;quit;
options ps=65 ls=171;

*/ /**************************************************************************************************************************/
*/ /*
*/ /*  TRANFORM BACK (SAME AS ORIGINAL PLOT
*/ /*
*/ /*
*/ /*       0        110       220       330       440       550       660       770       880
*/ /*       --+---------+---------+---------+---------+---------+---------+---------+---------+--
*/ /*    YT |                                                                                   |
*/ /*       |                                                                          ******   |
*/ /*  12.2 +   TRANSFORMATIONS BACK FROM UNIT SQUARE                        ***********        + 12.2
*/ /*       |                                                            ******                 |
*/ /*       |   xt=x *867;                                          ****                        |
*/ /*  11.2 +   yt = 2.22 + (12.49-2.22)*y ;                     ****                           + 11.2
*/ /*       |                                                    ***                            |
*/ /*       |                                                 ****                              |
*/ /*  10.2 +                                             ****                                  + 10.2
*/ /*       |                                         *****                                     |
*/ /*       |                                     *****                                         |
*/ /*   9.2 +                                   ***                                             +  9.2
*/ /*       |                                 ***                                               |
*/ /*       |                                **                                                 |
*/ /*   8.2 +                                * After fitting you want an estamte                +  8.2
*/ /*       |                                * in origial units                                 |
*/ /*       |                                *                                                  |
*/ /*   7.2 +                               *  Plug xmapped =657/867 = 0.7612456747  into       +  7.2
*/ /*       |                               *                                                   |
*/ /*       |                                  ymapped =sqrt(xmapped)*CDF('BETA',xmapped,a,b)   |
*/ /*   6.2 +                               *    + c*sin(1.57*mapped)                           +  6.2
*/ /*       |                               *                                                   |
*/ /*       |                                  y_back=miny+(maxy-miny)*ymapped;                 |
*/ /*   5.2 +                               *                                                   +  5.2
*/ /*       |                                  Check Model                                      |
*/ /*       |                               *                                                   |
*/ /*   4.2 +                               *  data _null_;                                     +  4.2
*/ /*       |                              **    x =657/867;                                    |
*/ /*       |                             **     y =2.22 + (sqrt(x)*CDF('BETA',x,307.3,500.4)   |
*/ /*   3.2 +                          ****         + 0.04218*sin(1.57*x))*(12.49-2.22);        +  3.2
*/ /*       |                     ******         put y=;                                        |
*/ /*       |          ************            run;quit;                                        |
*/ /*   2.2 + *********                                                                         +  2.2
*/ /*       |                                   Comparison Predicted vs Observed                |
*/ /*       |                                              X       Y                            |
*/ /*       |                                   Observed  657    11.76                          |
*/ /*       |                                  Predicted  657    11.56                          |
*/ /*       |                                                                                   |
*/ /*       |                                                                                   |
*/ /*       |                                                                                   |
*/ /*       |                                                                                   |
*/ /*       --+---------+---------+---------+---------+---------+---------+---------+---------+--
*/ /*         0        110       220       330       440       550       660       770       880
*/ /*                                           xt
*/ /*
*/ /**************************************************************************************************************************/


/*___                               _       _   _
|___ \   ___  _ __  ___   ___  ___ | |_   _| |_(_) ___  _ __
  __) | / _ \| `_ \/ __| / __|/ _ \| | | | | __| |/ _ \| `_ \
 / __/ | (_) | |_) \__ \ \__ \ (_) | | |_| | |_| | (_) | | | |
|_____| \___/| .__/|___/ |___/\___/|_|\__,_|\__|_|\___/|_| |_|
             |_|
*/


data  workx.check;

   a=8.8179917838  ;
   b=0.9573331691  ;
   c=341.0840722484;
   d=2.3689668227  ;
   e=1 ;
   set workx.rwant;
   py = A/(1 + B**(xs - C)) + D;
   row=_n_;
run;quit;

* join with ops solution;
options validvarname=v7;
proc sql;
  create
     table workx.addorg as
  select
     l.py as ops_y
    ,l.xs
    ,l.ys
    ,r.pt as Roger
  from                                           xt=x *867;
     workx.check as l, workx.back as r           yt =2.22 + (sqrt(x)*CDF('BETA',x,307.3,500.4)
  where                                                + 0.04218*sin(1.57*x))*(12.49-2.22);
     l.xs = r.xt
     and ys < 12.49
     and pt > 1
;quit;


options ls=100 ps=65;
proc plot data=workx.addorg;
 plot   ys*xs='.' ops_y*xs='o' roger*xs='r'/ overlay box
   haxis=0 to 900 by 100  vaxis=2.20 to 12.5 by 1;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*   POSTED SOLUTION ON STACKOVERFLOW COMPARED TO MY SOLUTION                                                             */
/*                                                                                                                        */
/*                                                                                                                        */
/*       --+---------+---------+---------+---------+---------+---------+---------+---------+---------+--                  */
/*       |                                                                                             |                  */
/*       |  OPS FIT ='o'                    ROGER = 'r'                   ORIGINAL DATA='.'            |                  */
/*       |                                                                                             |                  */
/*       |  a=8.8179917838                  A =  307.4                                                 |                  */
/*       |  b=0.9573331691                  B =  500.4                                                 |                  */
/*       |  c=341.0840722484                C = 0.0422                                                 |                  */
/*       |  d=2.3689668227                                                                             |                  */
/*       |  e=1 ;                                                                                      |                  */
/*       |                                                                                             |                  */
/*       |  py = A/(1 + B**(xs - C)) + D    Ymapped =sqrt(x)*CDF('BETA',x,a,b) + c*sin(1.57*x) [0,1]   |                  */
/*       |                                                                                             |                  */
/*       |                                  Final equation in original units                           |                  */
/*       |                                   yt =2.22 + (sqrt(x)*CDF('BETA',x,307.3,500.4) [0,900]     |                  */
/*       |                                          + 0.04218*sin(1.57*x))*(12.49-2.22);               |                  */
/*       |                                                                                             |                  */
/*       |                                                                                             |                  */
/*       | 0        100       200       300       400  XT   500       600       700       800       900|                  */
/*       --+---------+---------+---------+---------+---------+---------+---------+---------+---------+--                  */
/*    YT |                                                                                 rrrr ....   |                  */
/*       |                                                                               rr....        |                  */
/*  12.2 +  'o'=ops                                                                 ........           + 12.2             */
/*       |  'r'=rogers fit                                                     .....rr                 |                  */
/*       |  '.'=orginal data                                                ....rrr                    |                  */
/*       |                                                                ...rrr                       |                  */
/*       |                                                              ...rr                          |                  */
/*  11.2 +                                              oooooooooooooo...oooooooooooooooooooooooooo    + 11.2             */
/*       |                                           oooo           r..     stackoverflow fit          |                  */
/*       |                          'o'=ops fit-->  oo           rr..                                  |                  */
/*       |                                         oo          rr...                                   |                  */
/*       |                                        oo        rrr...                                     |                  */
/*  10.2 +                                        o       rr...                                        + 10.2             */
/*       |                                       oo     r....   <- '.'  original data                  |                  */
/*       |                                       o   rr...                                             |                  */
/*       |                                       o r....                                               |                  */
/*       |                                      or...                                                  |                  */
/*   9.2 +                                      o..                                                    +  9.2             */
/*       |                                     ...                                                     |                  */
/*       |                 'r'=rogers fit --> r.o                                                      |                  */
/*       |                                    .o                                                       |                  */
/*       |                                    .o                                                       |                  */
/*   8.2 +                                   .ro                                                       +  8.2             */
/*       |                                   .oo                                                       |                  */
/*       |                                   .o                                                        |                  */
/*       |                                   ro                                                        |                  */
/*       |                                   .o                                                        |                  */
/*   7.2 +                                   .o                                                        +  7.2             */
/*       |                                   o                                                         |                  */
/*       |                                  .o                                                         |                  */
/*       |                                   o                                                         |                  */
/*       |                                  ro                                                         |                  */
/*   6.2 +                                  .o                                                         +  6.2             */
/*       |                                  o                                                          |                  */
/*       |                                  .                                                          |                  */
/*       |                                  o                                                          |                  */
/*       |                                  o                                                          |                  */
/*   5.2 +                                 o.                                                          +  5.2             */
/*       |                                 or                                                          |                  */
/*       |                                 or                                                          |                  */
/*       |                                 o.                                                          |                  */
/*       |                                or.                                                          |                  */
/*   4.2 +                                or.                                                          +  4.2             */
/*       |                                or.                                                          |                  */
/*       |                               oo..  <- '.' Original data                                    |                  */
/*       |                               o .                                                           |                  */
/*       |                              o..r                                                           |                  */
/*   3.2 +                              ..r                                                            +  3.2             */
/*       |                           .... r                                                            |                  */
/*       |                          ..oo  r                                                            |                  */
/*       |                 .......  ooo  rr                                                            |                  */
/*       | ooooo............oooooo  rrrrrr                                                             |                  */
/*   2.2 + ......rrrrr                                                                                 +  2.2             */
/*       |                                                                                             |                  */
/*       --+---------+---------+---------+---------+---------+---------+---------+---------+---------+--                  */
/*         0        100       200       300       400       500       600       700       800       900                   */
/*                                                     XT                                                                 */
/**************************************************************************************************************************/
/*____                                        _
|___ /    ___ ___  _ __ ___  _ __   __ _ _ __(_)___  ___  _ __
  |_ \   / __/ _ \| `_ ` _ \| `_ \ / _` | `__| / __|/ _ \| `_ \
 ___) | | (_| (_) | | | | | | |_) | (_| | |  | \__ \ (_) | | | |
|____/   \___\___/|_| |_| |_| .__/ \__,_|_|  |_|___/\___/|_| |_|
                            |_|
*/
%utlfkil(d:/png/nlin_plots.png);

filename gout "d:/png/nlin.png";
ods html(id=add_dest) body=gout gpath="d:\png";
ods graphics on /reset=all imagename="nlin_plots";

title "pKa of an acid from titration";
proc sgplot data=workx.addorg;
label  ys="Original"
       ops_y="Stackoverflow"
;
 series x=xs y=roger    / lineattrs=(thickness=1pt color=blue);
 series x=xs y=ys    / lineattrs=(thickness=1pt  color=green);
 series x=xs y=ops_y / lineattrs=(thickness=1pt color=red);
 inset ( " " ="Stackoverflow fit "
         " " ="a=8.8179917838    "
         " " ="b=0.9573331691    "
         " " ="c=341.0840722484  "
         " " ="d=2.3689668227    "
         " " ="e=1 ;             "
         " " ="y = A/(1 + B**(xs - C)) + D"
 ) / position=topleft textattrs=( Size=9pt color=red) valuealign=left;
 inset ( " " ="Rogers fit "
         " " ="a =  307.4 "
         " " ="b =  500.4 "
         " " ="c = 0.0422 "
         " " ="e=1 ;      "
         " " ="Ymapped =sqrt(x)*CDF('BETA',x,a,b) + c*sin(1.57*x)"
 ) / position=bottomright textattrs=( Size=9pt color=blue) valuealign=left;
 /*--- inset ( /*--- " " ="Original Data "
 ) / position=right textattrs=( Size=9pt color=green) valuealign=left; ---*/
 inset ( " " ="Rogers fit "
         " " ="a =  307.4 "
         " " ="b =  500.4 "
         " " ="c = 0.0422 "
         " " ="e=1 ;      "
         " " ="Ymapped =sqrt(x)*CDF('BETA',x,a,b) + c*sin(1.57*x)"
 ) / position=bottomright textattrs=( Size=9pt color=blue) valuealign=left;

keylegend / location = inside position=topright;
run;quit;

filename gout clear;
ods graphics off;
ods _all_ close;
ods listing;

/*  _               _              _                           _ _ _          _
| || |     ___ __ _| | ___   _ __ | | ____ _   ___  __ _ _   _| (_) |__  _ __(_)_   _ _ __ ___
| || |_   / __/ _` | |/ __| | `_ \| |/ / _` | / _ \/ _` | | | | | | `_ \| `__| | | | | `_ ` _ \
|__   _| | (_| (_| | | (__  | |_) |   < (_| ||  __/ (_| | |_| | | | |_) | |  | | |_| | | | | | |
   |_|    \___\__,_|_|\___| | .__/|_|\_\__,_| \___|\__, |\__,_|_|_|_.__/|_|  |_|\__,_|_| |_| |_|
                            |_|                       |_|
*/

data tst;
  retain slope 0;
  do x=.001 to 1 by .001;
     xt=x *867;
     yt =2.22 + (sqrt(x)*CDF('BETA',x,307.3,500.4)  + 0.04218*sin(1.57*x))*(12.49-2.22);
     if dif(yt)/dif(xt) > slope then do; xs=xt; ys=yt;slope=dif(yt)/dif(xt);end;
     output;
  end;
  put xs= ys= slope=;
  stop;
run;quit;

/*--- xs=330.327 ys=5.7275695428 slope=0.1760288741 ---*/

options ls=64 ps=32;
proc plot data=tst;
 plot   yt*xt='*' / overlay box
   haxis=0 to 900 by 100  vaxis=2.20 to 12.5 by 1 href=330 vref=5.73;
run;quit;


              Plot of yt*xt.  Symbol used is '*'.

  yt --+-----+-----+-----+-----+-----+-----+-----+-----+--
12.2 +                     |                      *****  +
     |                     |                 ******      |
11.2 +                     |             *****           +
     |                     |          ****               |
10.2 +                     |      *****                  +
     |                     |   ****                      |
 9.2 +                     | ***                         +
     |                     |**                           |
 8.2 +                     |*                            +
     |                     **                            |
 7.2 +                     *                             +
     |                     *                             |
 6.2 + equlibrium pKa      *  (330,5.72)                 +
     |---------------------+-----------------------------|
 5.2 +                     *                             +
     |                    **                             |
 4.2 +                    *|                             +
     |                    *|                             |
 3.2 +                    *|                             +
     |                  ***|                             |
 2.2 + ******************  |                             +
     --+-----+-----+-----+-----+-----+-----+-----+-----+--
       0    100   200   300   400   500   600   700   800

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/

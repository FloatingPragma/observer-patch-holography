/* Derived adapter; original attribution retained below.
 * Corrections follow Martin PRD99 (2019) 033007 Eqs. (1.15), (2.16).
 * Original pinned SMDR source is never modified. */
/* The direct source for all results in this file: SPM 1812.04100,
   which is based in large part on earlier sources:

   Four-loop decoupling for alphaS, Schroeder and Steinhauser 0512058 and
   Chetyrkin, Kuhn, Sturm 0512060. Two-loop electroweak part from
   Bednyakov 1410.7603.

   Decoupling of alpha was given in eqs (6) and (7) of 9212285 Fanchiotti,
   Kniehl and Sirlin, and eq 29 of Kniehl, Pikelner and Veretin 1503.02138.

   QCD contribution to fermion mass decoupling at 2 loops and 3 loops
   can be found in eq. (27) of 0004189, based on eq. (20) of
   Chetyrkin, Kniehl, and Steinhauser 9708255.
*/

#include "smdr_internal.h"

/* --------------------------------------------------------------------- */
/* Combination is useful for finding fermion masses in the decoupled theory.
   qf = electric charge, I3f = weak isospin.
   Equation (2.20) of 1812.04100.
   For use internal to this file only.
*/
SMDR_REAL OPH_deltam1f (SMDR_REAL qf, SMDR_REAL I3f)
{
  return ((9.*g2 + 3.*gp2)/16. +
          qf * gp2 * (I3f + qf * (W/Z - 1.)) * (3.*SMDR_LOG(Z/Q2) - 2.5));
}

/* --------------------------------------------------------------------- */
/* Combination useful for matching relations for gauge couplings.
   Defined in equation (1.15) of 1812.04100.
   For use internal to this file only.
*/
SMDR_REAL OPH_Fthresh (SMDR_REAL x, SMDR_REAL y, SMDR_REAL qq)
{
  SMDR_REAL Ixxy, lnbarx, lnbary;

  Ixxy = TSIL_CREAL (TSIL_I2(x, x, y, qq));
  lnbarx = SMDR_LOG (x/qq);
  lnbary = SMDR_LOG (y/qq);
  return Ixxy + (x - y/2.) * lnbarx * lnbarx + y * lnbarx * lnbary +
         (4 * x - 2 * y) * lnbarx - 8 * x * lnbary +
         ((4*x - y) * (4*x - y)/(6 * x)) * (lnbary - lnbarx) -
         x/3. + (31./6.) * y - y * y/(3. * x);
}

/* ---------------------------------------------------------------------
   Simultaneously decouples the top quark and the Higgs, Z, and W bosons,
   at the current MSbar renormalization scale Q. The results are stored
   in the global variables:
   alphaS_53, alpha_53,	and m<fermion>_53,
   where <fermion> = b, c, s, u, d, tau, muon, electron.
   These are the MSbar QCD and EM couplings and the running masses in the
   QCD+QED effective theory that contains only the 5 quarks b,c,s,u,d and
   the charged leptons tau,mu,e, and the gluons and the photons.

   Currently the argument loopOrder is ignored, should fix that!!!
   --------------------------------------------------------------------- */
void OPH_Decouple_thZW (int loopOrder)
{
  /* char funcname[] = "SMDR_Decouple_thWZ"; */
  SMDR_REAL alphaS_nodec, alphaSoPi, alphaSoPi2, alphaSoPi3, alphaSoPi4;
  SMDR_REAL alpha_nodec, alphaoPi;
  SMDR_REAL lnbarT, lnbarW, lnbarZ, lnbarh, lnToW;
  SMDR_REAL lnbar2T, lnbar3T, lnbar4T, lnbar2W, lnbar2Z, lnbar2h;
  SMDR_REAL I0hW, I0hZ, I0TW, I0WZ, IhTT, IhTW, IhWW, IhZZ, ITTZ, ITWZ, IWWZ;
  SMDR_REAL eightWm5Z2, WmZ2, WmZ3, WmZ4, TmW, TmW2, TmW3, TmW4;
  SMDR_REAL theta2mbg34, theta2mbg32, theta2mbg30;
  SMDR_REAL theta2mfg34, theta2mfg32, theta2mfg30;
  SMDR_REAL Cf, I3f, Qf, Qf2, Qf3, Qf4;
  SMDR_REAL alpha_S_thresh_corr, alpha_thresh_corr;
  SMDR_REAL deltame1, deltame2, deltame;
  SMDR_REAL deltamu1, deltamu2g32, deltamu2g30, deltamu;
  SMDR_REAL deltamd1, deltamd2g32, deltamd2g30, deltamd;
  SMDR_REAL deltamb1, deltamb2g32, deltamb2g30, deltamb;
  SMDR_REAL deltamq2QCD, deltamq3QCD, deltamq4QCD;
  SMDR_REAL alphaS_dec_result, alpha_dec_result,
            mb_dec_result, mtau_dec_result,
            mc_dec_result, ms_dec_result, mmu_dec_result,
            mu_dec_result, md_dec_result, me_dec_result;

  /* QCD coupling at Q in full Standard Model with nothing decoupled. */
  alphaS_nodec = g32/(4 * PI);
  alphaSoPi = alphaS_nodec/PI;
  alphaSoPi2 = alphaSoPi * alphaSoPi;
  alphaSoPi3 = alphaSoPi2 * alphaSoPi;
  alphaSoPi4 = alphaSoPi3 * alphaSoPi;

  /* Fine structure constant at Q, in full Standard Model, nothing decoupled. */
  alpha_nodec = (g2 * gp2)/(4 * PI * g2pgp2);
  alphaoPi = alpha_nodec/PI;

  lnbarT = TVIL_LOG (T/Q2);
  lnbarW = TVIL_LOG (W/Q2);
  lnbarZ = TVIL_LOG (Z/Q2);
  lnbarh = TVIL_LOG (h/Q2);
  lnToW = lnbarT - lnbarW;

  lnbar2T = lnbarT * lnbarT;
  lnbar3T = lnbar2T * lnbarT;
  lnbar4T = lnbar3T * lnbarT;
  lnbar2W = lnbarW * lnbarW;
  lnbar2Z = lnbarZ * lnbarZ;
  lnbar2h = lnbarh * lnbarh;

  I0TW = TSIL_I2(0, T, W, Q2);
  I0hW = TSIL_I2(0, h, W, Q2);
  I0hZ = TSIL_I2(0, h, Z, Q2);
  I0WZ = TSIL_I2(0, W, Z, Q2);
  IhTT = TSIL_I2(h, T, T, Q2);
  IhTW = TSIL_I2(h, T, W, Q2);
  IhWW = TSIL_I2(h, W, W, Q2);
  IhZZ = TSIL_I2(h, Z, Z, Q2);
  ITTZ = TSIL_I2(T, T, Z, Q2);
  ITWZ = TSIL_I2(T, W, Z, Q2);
  IWWZ = TSIL_I2(W, W, Z, Q2);

  TmW = T-W;
  TmW2 = TmW * TmW;
  TmW3 = TmW2 * TmW;
  TmW4 = TmW3 * TmW;
  WmZ2 = (W-Z) * (W-Z);
  WmZ3 = WmZ2 * (W-Z);
  WmZ4 = WmZ3 * (W-Z);
  eightWm5Z2 = (8*W - 5*Z) * (8*W - 5*Z);

  /* Equation (2.20) of 1812.04100 */
  deltame1 = ONELOOPFACTOR * OPH_deltam1f (-1, -0.5);
  deltamu1 = ONELOOPFACTOR * OPH_deltam1f (2./3., 0.5);
  deltamd1 = ONELOOPFACTOR * OPH_deltam1f (-1./3., -0.5);

  /* Equation (2.21) of 1812.04100 */
  deltamb1 = deltamd1 + ONELOOPFACTOR * (3./4.) * yt2 * (
             W2 * lnToW/TmW2  - W/TmW + 5./6. - lnbarT);

  /* Equation (2.26) of 1812.04100. */
  deltamq2QCD = TWOLOOPFACTOR * g34 * (4./3.) * (
                89./36. + (5./3.) * lnbarT + lnbar2T);

  /* Equation (2.27) of 1812.04100. */
  deltamq3QCD = THREELOOPFACTOR * g36 * (
                126.16094650019454
                + 111.0479731067833 * lnbarT
                + (700./27.) * lnbar2T
                - (152./27.) * lnbar3T);

  /* Equation (2.28) of 1812.04100. */
  deltamq4QCD = FOURLOOPFACTOR * g38 * (
                236.908052
                + 452.388432 * lnbarT
                - 543.379386 * lnbar2T
                - (10984./81.) * lnbar3T
                + (830./27.) * lnbar4T);

  /* Equation (2.23) and ancillary file of 1812.04100. */
  #include "includes/theta2mb.c"
  deltamb2g32 = TWOLOOPFACTOR * theta2mbg32;
  deltamb2g30 = TWOLOOPFACTOR * theta2mbg30;
  deltamb = deltamb1 + deltamq2QCD + deltamq3QCD + deltamq4QCD +
            deltamb2g32 + deltamb2g30;

  /* Equation (2.26) and ancillary file of 1812.04100. */
  Cf = 4./3.; I3f = -0.5; Qf = -1./3.;
  #include "includes/theta2mf.c"
  deltamd2g32 = TWOLOOPFACTOR * theta2mfg32;
  deltamd2g30 = TWOLOOPFACTOR * theta2mfg30;
  deltamd = deltamd1 + deltamq2QCD + deltamq3QCD + deltamq4QCD +
            deltamd2g32 + deltamd2g30;

  Cf = 4./3.; I3f = 0.5; Qf = 2./3.;
  #include "includes/theta2mf.c"
  deltamu2g32 = TWOLOOPFACTOR * theta2mfg32;
  deltamu2g30 = TWOLOOPFACTOR * theta2mfg30;
  deltamu = deltamu1 + deltamq2QCD + deltamq3QCD + deltamq4QCD +
            deltamu2g32 + deltamu2g30;

  Cf = 0; I3f = -0.5; Qf = -1;
  #include "includes/theta2mf.c"
  deltame2 = TWOLOOPFACTOR * theta2mfg30;
  deltame = deltame1 + deltame2;

  mb_dec_result = (yb * v/SQRT2) * (1. + deltamb);
  mc_dec_result = (yc * v/SQRT2) * (1. + deltamu);
  ms_dec_result = (ys * v/SQRT2) * (1. + deltamd);
  md_dec_result = (yd * v/SQRT2) * (1. + deltamd);
  mu_dec_result = (yu * v/SQRT2) * (1. + deltamu);
  mtau_dec_result = (ytau * v/SQRT2) * (1. + deltame);
  mmu_dec_result = (ymu * v/SQRT2) * (1. + deltame);
  me_dec_result = (ye * v/SQRT2) * (1. + deltame);

  /* 4-loop QCD threshold correction for alphaS, for 6 quark -> 5 quark
     decoupling, from Schroeder and Steinhauser, 0512058.
     Three-loop and four-loop QCD parts are very small, but included anyway.
     Also see Chetyrkin, Kuhn, Sturm 0512060.
     The form used here is from 1812.04100 equations (2.15)-(2.18).
  */
  alpha_S_thresh_corr = 1.;

  alpha_S_thresh_corr += ONELOOPFACTOR * g32 * (2./3.) * lnbarT;

  alpha_S_thresh_corr += TWOLOOPFACTOR * g34 * ((4./9.) * lnbar2T
                     + (22./3.) * lnbarT + 22./9.);

  alpha_S_thresh_corr += THREELOOPFACTOR * g36 * ((8./27.) * lnbar3T
                     - 3. * lnbar2T + (620./9.) * lnbarT + 35.123151);

  alpha_S_thresh_corr += FOURLOOPFACTOR * g38 * ((16./81.) * lnbar4T +
                     (4706./81.) * lnbar3T - (1231./27.) * lnbar2T +
                     245.856958 * lnbarT - 109.765121);



  /* Next two contributions are from eq. (2.15) of 1812.04100,
     and equivalent to the earlier paper Bednyakov 1410.7603, eqs. (16)-(17).
     Quite small.
  */
  alpha_S_thresh_corr += TWOLOOPFACTOR * g32 * yt2 * (
     (2.*T*(h - T))/(h*(4.*T - h)*(4.*T - h)) * OPH_Fthresh (T, h, Q2)
     - (2./3.) * (T/h) * (1 + lnbarh - lnbarT)
     + (2./3.) * (T/((4.*T-Z)*(4.*T-Z)*(4.*T-Z))) * (
       T * (80.*W/Z - 7. - 64.*W2/Z2) + 8.*Z - 40.*W + 32.*W2/Z) *
       OPH_Fthresh (T, Z, Q2)
     + ((80.*W/Z - 64.*W2/Z2) * (1. + lnbarZ - lnbarT) + 2. +
       3.*lnbarh - 7.*lnbarZ - 14.*lnbarT)/18.);

  alpha_S_thresh_corr += TWOLOOPFACTOR * g32 * g2 * (
      ((8./9.)*(W-Z)/Z) * lnbarT + 3. * lnbarW
      + ((25./18.) * Z/W - 13./9. + (14./9.) * W/Z) * lnbarZ
      + (T/TmW) * lnToW - 49./27. - W/(18.*Z) - (163./216.) * Z/W);


  alphaS_dec_result = alphaS_nodec * alpha_S_thresh_corr;

  /* Next few results from 1812.04100 equations (2.13) and (2.14) */
  alpha_thresh_corr = 1.;

  alpha_thresh_corr += (alphaoPi/4) * (2./3. - 7.*lnbarW + (16./9.)*lnbarT);

  alpha_thresh_corr += (alphaoPi * alphaSoPi/16.) *
                       (-(64./9.) * lnbarT - 208./27.);

  alpha_thresh_corr += (alphaoPi/4.) * yt2 * ONELOOPFACTOR * (
    ((16.*T*(h-T))/(3.*h*(4.*T-h)*(4.*T-h))) * OPH_Fthresh(T, h, Q2)
    - (16.*T)/(9.*h) * (1. + lnbarh - lnbarT)
    + (16.*T)/(9.*(4.*T-Z)*(4.*T-Z)*(4.*T-Z)) * (T * (80.*W/Z - 7. -
       64.*W2/Z2) + 8.*Z - 40.*W + 32.*W2/Z) * OPH_Fthresh(T, Z, Q2)
    + 4.*( (80.*W/Z - 64.*W2/Z2) * (1. + lnbarZ - lnbarT) +
      2. + 3.*lnbarh - 7.*lnbarZ - 14.*lnbarT)/27. + 22.*lnbarT - 43./4.);

  alpha_thresh_corr += (alphaoPi/4.) * g2 * ONELOOPFACTOR * (
    (3. * W * (3. * h2 - 12. * h * W + 4. * W2)/
        (h * (4.*W - h) * (4.*W - h) * (4.*W - h))) * OPH_Fthresh(W, h, Q2)
    + (W/h - 2.) * lnbarh
    + (9. * W * (4. * W2 - 4. * W * Z + 3. * Z2)/
        (Z2 * (4.*W - Z) * (4.*W - Z))) * OPH_Fthresh(W, Z, Q2)
    + ((661. * Z)/(108. * W) - 491./27. + (319. * W)/(27. * Z) +
       (12. * W2)/Z2) * lnbarZ
    + (20./3. + (37. * W)/(3. * Z) - (12. * W2)/Z2 - W/h) * lnbarW
    + ((5. * T)/(3. * TmW)) * lnToW
    + 31./81. - (3. * h)/(4. * W) + W/h + (12. * W2/Z2) - (799. * W)/(27. * Z)
    - (1057. * Z)/(324. * W)) +
    (alphaoPi/4.) * (alphaoPi/4.) * (49 * lnbar2W -
      (224./9.) * lnbarT * lnbarW + (256./81.) * lnbar2T);

  alpha_dec_result = alpha_nodec * alpha_thresh_corr;

  Q_53 = Q;
  alphaS_53 = alphaS_dec_result;
  alpha_53  = alpha_dec_result;
  mb_53     = mb_dec_result;
  mc_53     = mc_dec_result;
  ms_53     = ms_dec_result;
  mu_53     = mu_dec_result;
  md_53     = md_dec_result;
  mtau_53   = mtau_dec_result;
  mmuon_53  = mmu_dec_result;
  melectron_53 = me_dec_result;

  return;
}

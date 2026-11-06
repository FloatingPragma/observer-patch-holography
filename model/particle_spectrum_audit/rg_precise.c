/* Precise replay of the declared SMDR RG boundary.
 * usage: rg_precise Qfinal loops inputfile [yt_override lambda_override]
 * Prints the final lambda, beta_lambda, yt, beta_yt, gauge couplings and v with 18 digits,
 * so that the criticality conditions lambda = beta_lambda = 0 can be checked beyond the
 * six printed digits of calc_RGrun. Links against the pinned SMDR library.
 */
#include "smdr_internal.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  char *varList[] = {"SMDR_Q_in","SMDR_lambda_in","SMDR_g3_in","SMDR_g_in","SMDR_gp_in",
    "SMDR_yt_in","SMDR_yb_in","SMDR_yc_in","SMDR_ys_in","SMDR_yu_in","SMDR_yd_in",
    "SMDR_ytau_in","SMDR_ymu_in","SMDR_ye_in","SMDR_v_in","SMDR_m2_in","SMDR_Lambda_in"};
  if (argc < 4) { fprintf(stderr, "usage: rg_precise Qfinal loops inputfile [yt lambda]\n"); return 2; }
  SMDR_REAL Qf = strtold(argv[1], 0);
  int loops = atoi(argv[2]);
  SMDR_Read_Values(argv[3], 17, varList);
  if (argc >= 6) {
    SMDR_yt_in = strtold(argv[4], 0);
    SMDR_lambda_in = strtold(argv[5], 0);
    SMDR_m2_in = -SMDR_lambda_in*SMDR_v_in*SMDR_v_in;
  }
  SMDR_RGeval_SM(Qf, loops);
  SMDR_Betas(loops);
  printf("{\"Q\":%.18Lg,\"lambda\":%.18Lg,\"beta_lambda\":%.18Lg,\"yt\":%.18Lg,\"beta_yt\":%.18Lg,"
         "\"g3\":%.18Lg,\"g2\":%.18Lg,\"gY\":%.18Lg,\"v\":%.18Lg,\"beta_v\":%.18Lg,\"loops\":%d}\n",
         Q, k, beta_k, yt, beta_yt, g3, g, gp, v, beta_v, loops);
  return 0;
}

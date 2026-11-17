# Preregistered source conservation/current join

Before any execution, bind this document, DESIGN.md and headline.py in PREREG.json.

1. Use the unchanged full_domain.branch_action for C and E and both coins.
   Enumerate every Fock input on two sites with two flavors and on four sites
   with two flavors, using seams (0,1) and (0,3), respectively. Retain all32
   coin/raw slots, including zeros and collisions. All arithmetic in the exact
   audit uses integer matrices (the Y generator is i times a real integer matrix).
2. Prove complete effects and unitality by exact matrix sums. Check every
   branch commutes with each doubled total spin generator and total number.
3. Check the exact conditional local continuity equation for every branch:
   V_h^T T_a V_h-T_a E_h=-j_h E_h(T_a+T_b), with the opposite sign at b;
   E_h=V_h^T V_h, K_h=V_h/sqrt(2), T=2S. Check every spectator spin separately.
4. Verify the su(2) algebra and on N=1 the exact Casimir S^2=3/4 and
   component spectrum +/-1/2. Check the polarized spatial trace is stationary
   under every seam in the twelve-port one-particle source. Its mean charge
   current must be zero. Do not interpret zero current as a magnetic moment.
5. Execute the existing native complete source dilation for both C/E on all16
   Fock inputs and all512 output addresses. Bound Frobenius map error by2e-7
   and conservation intertwiner error by2e-6. These are binary64 diagnostics,
   not outward-rounded or hardware certificates. Do not discard padding.
6. Execute global flavor swap using source-compiled fermionic swap words.
   Its ideal complete law must preserve every native coin/raw effect and
   current. The Sz conservation defect must be nonzero. Check its native
   matrix error below2e-7 and the sign-reversing intertwiner below2e-6.
7. Source-copy all four occupation bits to four fresh pointer bits after the
   repair and after the changed-law operation. Retain complete joint laws on
   all8192 field/coin/raw/pointer addresses for all16 initial basis states.
   Planted single-up and single-down must yield Sz=+/-1/2 with wrong-sign
   probability below2e-12; the changed-law control must reverse those laws.
   Raw normalization error must remain below2e-7. No conditional renormalization.

All gates fixed before execution. The new result is not a physical electron
identification, source-selected spin representation, Maxwell field moment,
magnetic gyromagnetic prediction or full completion of either row5 or6.

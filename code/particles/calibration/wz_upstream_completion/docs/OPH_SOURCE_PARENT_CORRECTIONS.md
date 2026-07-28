# OPH source-parent corrections required before an OPH-native W/Z action claim

These defects are upstream of the perturbative Standard-Model calculation. They do not block an **imported-Minkowski / imported-SM validation lane**, but they block an `OPH_NATIVE_PHYSICAL` claim.

## 1. Separate event position from observer-frame fiber

The correctly typed geometry is

```text
EventBase              M
Event                   x in M
FutureFrameFiberAt(x)   H3_frame(x) = {u in T_x M : g(u,u)=-1, u future}
ObserverFrame           u in H3_frame(x)
RestSpace               E_(x,u)=u^perp subset T_x M
```

The local Standard-Model action is integrated over `M`. A field is a section of a bundle over `M`. The hyperboloid `H3_frame(x)` parametrizes timelike frames at one event; it is not the event-position space.

### Mandatory mutation

Any packet containing

```text
record_position X_i(t) in H3_frame
```

without an independently constructed event-base map must fail `SPACETIME_PARENT`.

## 2. Mixed-GNS transport is not a MaxEnt corollary

A finite MaxEnt/refinement clause may say that coarse states remain in a finite exponential family. It does not produce:

- isometries/intertwiners between different GNS Hilbert spaces;
- compatible cyclic and separating vectors;
- compact-time convergence of modular automorphisms;
- inverse and group-law control; or
- support-covariant limiting modular flow.

Add an independent `MGNS-1` receipt with comparison maps `J_sr` and verify

```text
J_sr* J_sr = I
state/vector compatibility
compact-time Cauchy control for modular unitaries and inverses
identity and group law
support covariance
summable/cofinal residual modulus
```

Do not cite the MaxEnt axiom as if it contained these clauses.

## 3. Event-manifold parent

Before fields can be local functions/sections, the event branch must provide:

- a Hausdorff four-manifold or a certified local chart packet;
- local surjectivity/homeomorphism, not merely density of record images;
- overlap cocycles;
- an oriented time-oriented Lorentz metric; and
- a spin structure.

A complete embedded metric subspace of `R^4` need not be open. Use a local-degree/surjectivity certificate or state the manifold as an explicit branch input.

## 4. Lorentz-cone parent

A projectivized boundary conformal to `S^2` does not by itself force the cone to be quadratic. Add `CONE-QUADRIC-1`:

```text
nondegenerate symmetric bilinear form g_x
boundary equals {v != 0 : g_x(v,v)=0}
inertia exactly (1,3) or (3,1)
continuous/refinement-compatible orientation
```

## 5. Clock parent

Relative gamma

```math
-g(u,u')
```

is a local relative-velocity observable, not a universal transition map between two proper-time histories. A clock receipt needs worldlines, event correspondence/synchronization, and

```math
d\tau^2=-g_{ab}dx^a dx^b
```

with calibration and uncertainty. This geometric clock parent is separate from the source calculation of the operational reference transition used to attach GeV units.

## 6. Consequence for the W/Z program

Use two lanes:

1. `EXTERNAL_SM_EFT_VALIDATION`: explicitly imported Minkowski chart and complete imported SM packet. This validates FJ, BRST, Nielsen, analytic continuation and pole extraction.
2. `OPH_NATIVE_*`: additionally requires the corrected event/manifold/cone/MGNS/clock parent receipts plus the source action and matching packet.

Never let success of lane 1 mark an OPH geometry or source-law gate true.

# Parallel string stress and fixed-volume ball areas

A source with positive energy and longitudinal tension has a curved transverse geometry in this version of partition repair. The released traveller and home bits specify the ordinary and correlation sheet densities through the traveller count and their same-service coincidence count. A transverse rod survey calibrates the angle defect per unit ordinary sheet energy. At the same gain the second density predicts the signed area change of ten moving observers through one event. The model supplies a string action, a transport rule and proper rods and clocks. The calculations below concern that extension. [Action, transport law and independent survey](geometry.py) · [Finite-radius measurements](geometry.json)

## What the material and geometric readers measure

Use units in which the rods, clocks and sheet tension are one. Let the two densities $o,c$ each lie in $[0,1]$. A frozen lottery on the archived services registers the homogeneous proper transverse number density $\rho=o+c$. The geometric rule assigns rotation $2\pi$ per unit integrated sheet number. In a positive local orthonormal $(x,y)$ plane its curvature two-form is $2\pi\rho\,dA$. Parallel sheets run along the time and $z$ directions. The extension gives observers proper radial rods, angular bearings, transport of a transverse direction and Minkowski clocks along each sheet. On a transverse disk about the common crossing the rotation rule integrates to
\[
 ds^2=-dt^2+dz^2+dr^2+S(r)^2 d\phi^2,\qquad
 S''=-K S,\quad S(0)=0,\quad S'(0)=1,\quad K=2\pi(o+c).
\]
Indeed the oriented rotation around a circle is $2\pi(1-S'(r))$, and its enclosed proper area is $2\pi\int_0^r S(s)ds$. A rod of length $r$ gives $S(r)=\sin(\sqrt K r)/\sqrt K$ or $r$ for $K=0$. A transported bearing gives the rotation, a circumference gives $2\pi S(r)$, and radial rings give the area. The input $K$ uses the registered densities alone. Neither the stress inferred by an observer nor a ball entropy or an Einstein equation enters the radial transport rule.

A disk of rod radius $1/8$ on the ordinary-only preparation gives rotation $0.10252832$ radians and area $0.04895367$. Dividing these readings by each other and by the independently registered ordinary energy density $1/3$ gives $\widehat\kappa=6.28318530718$. Freeze this value before changing the second record. This is a calibration between independent instruments in a designed geometry. The density is a locally uniform continuum material introduced by the extension. Concentrating a full $2\pi$ defect at an isolated thin point, deriving a smooth density limit from finitely many such points, and assembling an entire compact transverse surface are different constructions.

## Action, conservation and tensor components

Each sheet has action $-\int\sqrt{-\det\gamma}\,dt\,dz$, where $\gamma$ is the pullback of the spacetime metric to that sheet. Sum these actions against a conserved transverse material-label measure $d\nu$, with background representation $d\nu=\rho\,dA$. Hold $d\nu$ fixed when varying the metric. The longitudinal force measured by slowly lengthening a straight sheet is its tension; the energy per proper longitudinal length is the same unit. Both can be determined without a moving-ball area. Varying $g_{ab}$ gives $\delta I=\tfrac12\int T^{ab}\delta g_{ab}\,d\mathrm{vol}_4$. In the local orthonormal frame the covariant stress components are
\[
 T_{ab}=\rho(U_a U_b-Z_a Z_b)=\operatorname{diag}(\rho,0,0,-\rho).
\]
Here $U$ is the unit timelike sheet direction and $Z$ the unit tangent along $z$. Each family has its own conserved label measure. The sheets are totally geodesic, $U,Z$ are parallel in the product, and both directional derivatives of $\rho$ vanish. These facts give $\nabla^aT_{ab}=0$ for each family independently. An observer moving with unit velocity $u$ measures $T_{uu}=\rho[(u^0)^2-(u^z)^2]=\rho b^2$, where $b^2=1+(u^x)^2+(u^y)^2$.

Compute curvature separately using the rod metric. The surface Ricci tensor is $K h_{xy}$. The flat directions have zero Ricci tensor, and scalar curvature is $2K$. Subtracting half the scalar times the metric therefore gives $G_{ab}=\operatorname{diag}(K,0,0,-K)$ and $G_{uu}=Kb^2$. Thus the common ordinary calibration $\kappa=2\pi$ gives $G_{ab}=\kappa T_{ab}$ on this nonvacuum homogeneous sector with flat empty reference and $\Lambda=0$. The same calibration means $\eta=2\pi/\kappa=1$ in these units. A factor $1/2$ or zero on the correlation record in the transport rule leaves the sheet action intact and fails the held comparison.

## Independent proper balls in ten frames

The common event is the origin. In the declared product chart the observers can realize their spacetime directions by rods, clocks and the local Lorentz metric. For each of the ten registered velocities, take $\exp_p(u^\perp)$ and stop its proper geodesic ball at a radius $r$ that encloses $4\pi\ell^3/3$. The radial Gauss lemma holds on this slice, and its extrinsic curvature at the central event vanishes. The three-dimensional Gauss equation gives $R^{(3)}=2G_{uu}=2Kb^2$. In the table the direct integral of the boundary metric is evaluated at its separately solved volume stop. It uses neither the action readings nor the Gauss prediction. The paired $z$ boosts see the rest result. Of the mixed boosts the $xy$ observer has the largest $b^2=9$. [Ten observer directions and inverse](../../../../code/matched_einstein_limit/frame.py) · [Existing exponential-map quadrature](../../../../code/matched_einstein_limit/finite_entropy_region/refinement.py)

The entries give $(A_{\rm held}-A_{\rm ordinary})/\ell^4$ at $\ell=1/256$. The forecast is $-(4\pi/15)\widehat\kappa\,(T^{\rm held}-T^{\rm ordinary})_{uu}$. On the held lottery the traveller and home single probabilities are fixed and the coincidence density increases from zero to $8/33$. There are no fitted parameters in these held transfers.

| Observer | $b^2$ | Stress forecast | Proper area | Half-strength angular rule | Unread correlation |
|---|---:|---:|---:|---:|---:|
| Rest, $\pm z$ | 1 | −1.276070 | −1.276074 | −0.638037 | 0 |
| $\pm x,\pm y$ | 25/9 | −3.544639 | −3.544703 | −1.772347 | 0 |
| $xz,yz$ | 5 | −6.380350 | −6.380586 | −3.190277 | 0 |
| $xy$ | 9 | −11.484631 | −11.485444 | −5.742668 | 0 |

These are deterministic model integrals in proper rod units. The calculation covers all ten directions at four scales and archives both actual radii, both areas, the stress-weighted proper volumes and both controls. It tests source conservation and the constitutive coupling on a homogeneous family of parallel sheets. The supplied angle-defect rule encodes the local strength of gravity in the sector. An independent derivation of that rule or of the string action from freely evolved Interval Echo matter, an observational galaxy or lens prediction, and general gravitational dynamics lie outside this result.

## Uniform finite-radius errors

Write $H=K b^2$, $\sigma=H\ell^2$. Our domain $K\le4\pi$, $1\le b^2\le9$, $\ell\le1/32$ gives $0\le\sigma<1/8$. For the rest-slice unit direction $(\mu,\phi)$ of the quadrature, put $f=b^2\cos^2\phi+\sin^2\phi$ and $q^2=(1-\mu^2)f$. The angular and volume Jacobian is
\[
 J=\sqrt{1+d},\qquad d={b^2\over f}\left[\operatorname{sinc}^2(\sqrt K r q)-1\right].
\]
It comes from the single transverse angular Jacobi field in the curved factor. The direction parallel to the ray has unit stretch. For $r\le2\ell$, let $\alpha=Hr^2\le1/2$, $w=\alpha(1-\mu^2)$. The cosine series for $\operatorname{sinc}^2 x=(1-\cos(2x))/(2x^2)$ gives $-w/3\le d\le0$ and $|d+w/3|\le(2/45)\alpha^2$. The identity $\sqrt{1+d}-1-d/2=-d^2/[2(1+\sqrt{1+d})^2]$ then gives $|J-1+w/6|\le\alpha^2/12$. Its angular mean has leading term $1-Hr^2/9$ since the average of $1-\mu^2$ is $2/3$.

In units of the target volume, $F(t)=V(\ell t)/(4\pi\ell^3/3)$ obeys
\[
 |F(t)-t^3+\sigma t^5/15|\le\sigma^2t^7/28,\qquad 0\le t\le2.
\]
Also $F(1)\le1$, $F(2)>1$ and $F'(t)=3t^2\langle J\rangle>2$ on $[1,2]$. For $t_0=1+\sigma/45$ use $t_0<1.01$ to bound $|F(t_0)-1|/\sigma^2$ by $(3+1/360)/45^2+5(1.01)^4/(15\cdot45)+(1.01)^7/28<0.05$. The exact stop $t=r/\ell$ is therefore within $0.025\sigma^2$ of $t_0$, and $0\le t-1\le\sigma/30$. This includes the flat case by exact evaluation. Inserting the radius in $A/(4\pi\ell^2)=t^2-\sigma t^4/9+E$, $|E|\le\sigma^2t^6/12$, gives
\[
 \left|A-4\pi\ell^2+{4\pi\over15}H\ell^4\right|\le 4\pi H^2\ell^6 .
\]
For example, the normalized errors in $t^2$, $\sigma t^4/9$ and $E$ after their terms through order $\sigma$ are below $0.052\sigma^2,0.016\sigma^2,0.089\sigma^2$ respectively. The bound is uniform on the full density range and all ten frames. The errors add for contrasts between two nonempty sectors.

The radial stress weight on the real slice is $W=\int_{B(r)}(r^2-s^2)/(2r)\,dV=\int_0^r(s/r)V(s)ds$. The coarse inequality $|J-1|\le Hs^2/3$ gives $|W-4\pi r^4/15|\le4\pi H r^6/63$. At the matched stop, $t<1.01$ and $|t-1|\le\sigma/30$ give
\[
 |W-4\pi\ell^4/15|\le (4\pi/20)H\ell^6.
\]
Thus an independently registered uniform sheet stress has weighted action readout $2\pi T_{uu}W$ with a controlled sixth-order geometric bias. Identifying this readout with the entropy of particular finite words requires their own state and dynamics; that identity is not assumed in this geometry calculation. [Checks of the derivatives, volume stops, action and second weighted integration](test_geometry.py)

~~~sh
python -m model.candidate_laws.interval_echo_partition_repair.published_strings.geometry
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/published_strings/test_geometry.py
~~~

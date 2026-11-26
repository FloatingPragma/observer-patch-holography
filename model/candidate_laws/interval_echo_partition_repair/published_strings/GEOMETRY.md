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

## A finite survey using neighboring rods

The same held comparison can be read without calling the smooth area or volume integral. Along a transverse ray of total proper length $s$, place $N=\lceil 1/\ell\rceil$ adjacent radial rods of length $h=s/N$. Start the signed transverse separation at $J_0=0$, $J_1=h$. The material angle rule gives $J_{m+1}=2J_m-J_{m-1}-Kh^2J_m$. Two neighboring longitudinal differences and the locally enrolled angular defect determine the next pointer. This is an additional finite transport rule for the published sheets. Its geometry input is the density-dependent defect strength; supplying a continuum ball area to the procedure would make the test redundant.

For $z=Ks^2$ the unrounded final ratio is the polynomial
\[
 {J_N\over s}=\sum_{k=0}^{N-1}{(-z)^k\over(2k+1)!}
         \prod_{j=1}^k(1-j^2/N^2).
\]
Indeed the coefficient before division by $N$ in $J_N/h$ is $\binom{N+k}{2k+1}$ for powers of $-Kh^2$, by the recurrence and Pascal's rule. The alternating terms decrease for $0\le z\le1/2$. Thus $|J_N/s-1+(1-N^{-2})z/6|\le z^2/120$. The reader rounds each recurrence increment in units of $N^{-8}$ of the total radial length. An additive error at one step has at most $m$ times its amplitude after $m$ steps, since the recurrence impulse response is $\sin(m\theta)/\sin\theta$ with $2\cos\theta=2-z/N^2$. The bound follows by summing $m$ unit complex phases and includes $z=0$ by continuity. Accumulated rod error is at most $N^{-6}$. The axis-square-root pointer rounds down by less than $N^{-8}$. The smooth density stays above $\sqrt{5/6}$, so that square root is uniformly Lipschitz on the perturbed readings. For rational density registers the alternating Machin formula for $2\pi$ supplies a finite angle calibration with error below $\ell^4/16$.

Take the six positive and negative principal rays of the transverse projection with equal positive weights. This cubature integrates every quadratic on the unit sphere exactly. Along their three axes the finite angular densities are, respectively, the rod ratio for $z=Kb^2r^2$, the square root of $1+b^2[(J_N/s)^2-1]$ at $z=Kr^2$, and one. Their mean is $1-(1-N^{-2})Hr^2/9+O(H^2r^4+N^{-6})$. Seven rational radial nodes with positive weights integrate each power through degree seven exactly. Relative to the true integral at any $r\le2\ell$, the difference in mean angular density is $O(Hr^2/N^2+H^2r^4+N^{-6})$. Consequently the volume has absolute error $O(\ell^7)$ uniformly over the enrolled bounded $K,b^2$. Stop against the rational target $\ell^3/3$ with the common solid-angle factor removed. Bisection localizes the final radius to a bin at most $\ell^5$ wide. The true volume derivative on $[\ell/2,2\ell]$ is bounded below by a constant times $\ell^2$. Every erroneous side decision is therefore within $O(\ell^5)$ of the true matched radius, which bounds the final radius error even with rounded readings. At that stop the inferred area differs from the ideal proper area by $O(\ell^6)$ and the weighted volume by $O(\ell^8)$. These bounds apply through refinement with growing finite pointer precision.

On the archived $8/33$ coincidence increase the $xy$ area contrast per $\ell^4$ is $-11.536904$, $-11.486730$ and $-11.484776$ at $\ell=1/32,1/128,1/512$ respectively. The frozen stress forecast is $-11.484631$. All ten calibrated directions can be tested on the same radii. In reporting physical area the output multiplies the rational area register by the common solid-angle unit $4\pi$. [Finite rod and stopping registers](finite_rods.py) · [Polynomial, rounding and all-frame comparisons](test_finite_rods.py)

~~~sh
python -m model.candidate_laws.interval_echo_partition_repair.published_strings.geometry
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/published_strings/test_geometry.py
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/published_strings/test_finite_rods.py
~~~

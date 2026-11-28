# An all-distance envelope for the selected pair

Let Δⱼ be the *signed* drop at the outer edge rⱼ of each raw annulus, including the final drop to zero. Thus the measured radial density is the sum of the filled disks of radii rⱼ with coefficients Δⱼ. Under a change of distance by a factor d, the radii scale by d and the surface densities stay fixed. Apparent integrated luminosities and gas masses scale as d². At fixed angle, the ordinary squared speed, including the exterior molecular correction, scales as d. The line-of-sight orbital speed is held fixed when re-deprojecting by a new inclination. Source profiles stay at the published Leroy inclinations and shapes throughout the distance box. [Complete implementation](distance.py)

In an area-preserving coarse grid, a unit filled disk with edge r inside a coarse shell [a,b] has an inner plateau and an outer partial shell of density f=(r²−a²)/(b²−a²). The effective boundary radius is q=a+(b−a)f. Hence

$$
0\le r-q=\frac{(r-a)(b-r)}{a+b}\le\frac{w}{4},
$$

where w is the maximum bin width. Also a+b≥2r−w. The upper bound on r−q can therefore be tightened to w²/[4 max(w,2r−w)]. The coarse field preserves linear superposition. Its unsigned boundary variation is at most the sum of the individual unsigned disk variations, and is at least its signed radial jump sum. For every d∈[dₗ,dₕ] this gives the simultaneous enclosure

$$
2\pi\left[d_l\sum_j r_j\Delta_j-\sum_j (\Delta_j)_+\frac{w^2}{4\max(w,2d_l r_j-w)}\right]
\ \le S(d) \le\ 2\pi d_h\sum_j r_j|\Delta_j|.
$$

All raw surface densities of the pair stay below the threshold ceiling. The signed first moment is positive since it equals the radial integral of the input density. The bounds include partially filled last bins and changes to the number of coarse bins as distance varies. The numerical test checks both inequalities on 20 synthetic nonmonotone profiles over 61 distance factors, and verifies mass conservation. The velocity extrema over the boxes follow directly from monotonicity of the sine at the supplied inclinations and from the positive baryonic force. Dividing the positive extremes of Y and S supplies a conservative outer enclosure rather than an extremum inferred by sampling distance endpoints. [Checks](test_measure.py) · [One- and two-width boxes](distance.json)

To bound a missing smooth exterior tail, suppose its surface density is initially at most the interior edge density and is nonincreasing thereafter, with rΣ(r)→0. Integration by parts combines its internal jumps and its join at R with the removed terminal edge. The net additional variation is ΔS=2π∫ᴿ∞Σ(r)dr≤Mtail/R. For a general positive, constant-density annular component of width h, the two boundaries give perimeter per mass 2/h. By the triangle inequality, any sum of full annuli with width at least w adds at most 2M/w to the source. These inequalities concern the radial observed source and its exterior. Neither constrains fine azimuthal structure at fixed annular mass. [Inventory and mass requirements](distance.json)

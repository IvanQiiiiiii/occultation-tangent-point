# occultation-tangent-point
MATLAB implementation to compute COSMIC occultation tangent (impact) points from LEO–GPS ECEF positions using WGS-84 ellipsoid geometry.

This repository provides a MATLAB implementation to compute the **tangent point** of a COSMIC radio occultation event.  
It calculates the **closest point** of the line-of-sight between a LEO and a GNSS satellite to the Earth's center, then converts that point to **geodetic latitude, longitude, and altitude** using the WGS-84 ellipsoid.

---

## 🌍 Features

- Input: Two satellite positions in **ECEF coordinates (km)** — LEO and GPS.  
- Output:  
  - Tangent point radius (impact parameter `p`)  
  - Geodetic latitude & longitude of the tangent point  
  - Height above the reference ellipsoid (`h`)  
- Option to use:
  - **Radial height** (analytical method)  
  - **Geodetic height** (using MATLAB’s built-in `ecef2geodetic` for higher accuracy)  
- Supports `wgs84Ellipsoid('kilometer')` by default, or any user-provided reference ellipsoid.  
- Includes metadata (e.g., LOS direction, whether the tangent lies between the satellites, etc.).

---

## 🧮 Formula Summary

Given ECEF coordinates of LEO (`r₁`) and GPS (`r₂`):

$$
\mathbf{n} = \frac{\mathbf{r}_2 - \mathbf{r}_1}{\|\mathbf{r}_2 - \mathbf{r}_1\|}, \quad
\mathbf{r}_{\perp} = \mathbf{r}_1 - (\mathbf{r}_1 \cdot \mathbf{n})\mathbf{n}, \quad
p = \|\mathbf{r}_{\perp}\|
$$

$$
\lambda = \mathrm{atan2}(p_y, p_x), \quad
\phi = \arctan\!\left[\left(\frac{a}{b}\right)^{2} \frac{p_z}{\rho}\right], \quad
h = r\!\left(1 - \frac{ab}{\sqrt{a^{2}p_z^{2} + b^{2}\rho^{2}}}\right)
$$

or use MATLAB’s:
```matlab
[lat, lon, h] = ecef2geodetic(px, py, pz, wgs84Ellipsoid('kilometer'));
```

## 💻 Example
```matlab
% Example LEO & GPS ECEF positions (km)
x_LEO = 6.262497051476043e+03;
y_LEO = 9.048916507518953e+02;
z_LEO = 2.788150586327804e+03;
x_GPS = -1.035909043931053e+04;
y_GPS = 2.182968511286796e+04;
z_GPS = 8.190449811103248e+03;

[p, h, glat, glon, meta] = impact_TangentPoint(x_LEO, y_LEO, z_LEO, x_GPS, y_GPS, z_GPS, 'method','radial');

fprintf('Tangent point: %.2f km, lat %.2f°, lon %.2f°, height %.2f km\n', ...
        p, lat, lon, h);
```

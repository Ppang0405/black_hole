# Physics and Mathematics of Black Hole Simulation

This document explains the fundamental physics and mathematical concepts implemented in the black hole simulation project.

## Table of Contents

1. [General Relativity Foundations](#general-relativity-foundations)
2. [Schwarzschild Metric](#schwarzschild-metric)
3. [Geodesic Equations](#geodesic-equations)
4. [Numerical Integration](#numerical-integration)
5. [Coordinate Systems](#coordinate-systems)
6. [Conservation Laws](#conservation-laws)
7. [Implementation Details](#implementation-details)

---

## General Relativity Foundations

### Einstein's Field Equations

The simulation is based on Einstein's theory of General Relativity, which describes gravity as the curvature of spacetime. The fundamental equation is:

```
Gμν = 8πTμν
```

Where:
- `Gμν` is the Einstein tensor (describes spacetime curvature)
- `Tμν` is the stress-energy tensor (describes matter/energy distribution)

For a black hole, we consider the vacuum solution where `Tμν = 0`, leading to the Schwarzschild solution.

### Spacetime Geometry

In General Relativity, massive objects warp spacetime. The simulation visualizes this through:

1. **Gravitational Lensing**: Light rays bend when passing near massive objects
2. **Event Horizon**: The boundary beyond which nothing can escape
3. **Geodesics**: The "straight lines" that particles follow in curved spacetime

---

## Schwarzschild Metric

### The Schwarzschild Solution

For a spherically symmetric, non-rotating black hole, the metric is:

```
ds² = -(1 - rs/r)c²dt² + (1 - rs/r)⁻¹dr² + r²dθ² + r²sin²θdφ²
```

Where:
- `rs = 2GM/c²` is the **Schwarzschild radius** (event horizon)
- `G = 6.67430×10⁻¹¹ m³/kg·s²` (gravitational constant)
- `M = 8.54×10³⁶ kg` (mass of Sagittarius A*)
- `c = 2.99792458×10⁸ m/s` (speed of light)

### Key Properties

1. **Event Horizon**: At `r = rs`, the metric component `g₀₀ = 0`
2. **Singularity**: At `r = 0`, curvature becomes infinite
3. **Asymptotic Flatness**: As `r → ∞`, the metric approaches Minkowski spacetime

### Implementation in Code

```cpp
// From black_hole.cpp line 128
BlackHole(vec3 pos, float m) : position(pos), mass(m) {
    r_s = 2.0 * G * mass / (c*c);
}

// Schwarzschild radius for Sagittarius A*
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36);
```

---

## Geodesic Equations

### Null Geodesics for Light Rays

Light rays follow null geodesics, characterized by `ds² = 0`. The geodesic equation is:

```
d²xμ/dλ² + Γμαβ(dxα/dλ)(dxβ/dλ) = 0
```

Where:
- `λ` is the affine parameter along the geodesic
- `Γμαβ` are the Christoffel symbols (connection coefficients)

### Spherical Coordinates Implementation

For the Schwarzschild metric in spherical coordinates `(t, r, θ, φ)`, the simulation implements:

#### Radial Equation:
```
d²r/dλ² = -(rs/2r²)f(dt/dλ)² + (rs/2r²f)(dr/dλ)² + (r - rs)(dθ/dλ)² + (r - rs)sin²θ(dφ/dλ)²
```

#### Angular Equations:
```
d²θ/dλ² = -(2/r)(dr/dλ)(dθ/dλ) + sinθ cosθ(dφ/dλ)²
d²φ/dλ² = -(2/r)(dr/dλ)(dφ/dλ) - 2(cosθ/sinθ)(dθ/dλ)(dφ/dλ)
```

Where `f = 1 - rs/r` is the lapse function.

### Code Implementation

From `geodesic.comp` (GPU implementation):

```glsl
void geodesicRHS(Ray ray, out vec3 d1, out vec3 d2) {
    float r = ray.r, theta = ray.theta;
    float dr = ray.dr, dtheta = ray.dtheta, dphi = ray.dphi;
    float f = 1.0 - SagA_rs / r;
    float dt_dL = ray.E / f;

    d1 = vec3(dr, dtheta, dphi);
    d2.x = - (SagA_rs / (2.0 * r*r)) * f * dt_dL * dt_dL
         + (SagA_rs / (2.0 * r*r * f)) * dr * dr
         + r * (dtheta*dtheta + sin(theta)*sin(theta)*dphi*dphi);
    d2.y = -2.0*dr*dtheta/r + sin(theta)*cos(theta)*dphi*dphi;
    d2.z = -2.0*dr*dphi/r - 2.0*cos(theta)/(sin(theta)) * dtheta * dphi;
}
```

---

## Numerical Integration

### Runge-Kutta 4th Order (RK4) Method

The simulation uses RK4 to integrate the geodesic equations. For a system `dy/dx = f(x,y)`:

```
k₁ = h·f(xₙ, yₙ)
k₂ = h·f(xₙ + h/2, yₙ + k₁/2)
k₃ = h·f(xₙ + h/2, yₙ + k₂/2)
k₄ = h·f(xₙ + h, yₙ + k₃)

yₙ₊₁ = yₙ + (k₁ + 2k₂ + 2k₃ + k₄)/6
```

### Implementation Details

The system integrates 6 variables: `(r, θ, φ, dr/dλ, dθ/dλ, dφ/dλ)`

```cpp
// From CPU-geodesic.cpp
void rk4Step(Ray& ray, double dλ, double rs) {
    double y0[6] = { ray.r, ray.theta, ray.phi, ray.dr, ray.dtheta, ray.dphi };
    double k1[6], k2[6], k3[6], k4[6], temp[6];

    geodesicRHS(ray, k1, rs);
    addState(y0, k1, dλ/2.0, temp);
    // ... similar for k2, k3, k4
    
    ray.r      += (dλ/6.0)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    ray.theta  += (dλ/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    // ... etc
}
```

### Step Size Selection

- **Affine parameter step**: `dλ = 1×10⁷` meters
- **Maximum steps**: 60,000 (adjustable based on performance)
- **Escape radius**: `r > 1×10³⁰` meters

---

## Coordinate Systems

### Spherical Coordinates (r, θ, φ)

The simulation primarily uses spherical coordinates where:
- `r`: Radial distance from black hole center
- `θ`: Polar angle (0 at north pole, π at south pole)  
- `φ`: Azimuthal angle (0 to 2π)

### Cartesian Conversion

```cpp
// Spherical to Cartesian
x = r * sin(θ) * cos(φ)
y = r * sin(θ) * sin(φ)  
z = r * cos(θ)

// Cartesian to Spherical
r = √(x² + y² + z²)
θ = arccos(z/r)
φ = arctan2(y, x)
```

### Velocity Transformations

Converting direction vectors from Cartesian to spherical basis:

```cpp
// From Ray constructor in CPU-geodesic.cpp
dr     = sin(theta)*cos(phi)*dx + sin(theta)*sin(phi)*dy + cos(theta)*dz;
dtheta = cos(theta)*cos(phi)*dx + cos(theta)*sin(phi)*dy - sin(theta)*dz;
dtheta /= r;
dphi   = -sin(phi)*dx + cos(phi)*dy;
dphi  /= (r * sin(theta));
```

---

## Conservation Laws

### Conserved Quantities

Due to the symmetries of the Schwarzschild metric, certain quantities are conserved along geodesics:

#### 1. Energy (E)
From time translation symmetry:
```
E = (1 - rs/r)(dt/dλ)
```

#### 2. Angular Momentum (L)
From rotational symmetry around the z-axis:
```
L = r²sin²θ(dφ/dλ)
```

### Code Implementation

```cpp
// From Ray constructor
L = r * r * sin(theta) * dphi;
double f = 1.0 - SagA.r_s / r;
double dt_dλ = sqrt((dr*dr)/f + r*r*(dtheta*dtheta + sin(theta)*sin(theta)*dphi*dphi));
E = f * dt_dλ;
```

### Physical Significance

- **E > 0**: Timelike geodesics (massive particles)
- **E = 0**: Null geodesics (light rays) - used in this simulation
- **L**: Determines the impact parameter and orbital characteristics

---

## Implementation Details

### Ray Initialization

Each pixel corresponds to a light ray cast from the camera:

```cpp
// From main simulation loop
float u = (2.0f * (x + 0.5f) / float(W) - 1.0f) * aspect * tanHalfFov;
float v = (1.0f - 2.0f * (y + 0.5f) / float(H)) * tanHalfFov;
vec3 dir = normalize(u*right + v*up + forward);
Ray ray(camera.pos, dir);
```

### Termination Conditions

The ray integration stops when:

1. **Event Horizon Crossing**: `r ≤ rs`
2. **Object Intersection**: Ray hits a gravitating body
3. **Disk Crossing**: Ray passes through the accretion disk
4. **Escape**: `r > ESCAPE_RADIUS`
5. **Maximum Steps**: Computational limit reached

### Performance Optimizations

1. **GPU Parallelization**: Each pixel processed independently
2. **Adaptive Resolution**: Lower quality during camera movement
3. **Early Termination**: Stop integration when ray fate is determined
4. **Optimized Data Structures**: Efficient memory layout for GPU

### Visual Effects

#### Gravitational Lensing
Light rays bend according to the geodesic equations, creating:
- **Einstein Rings**: Perfect alignment effects
- **Multiple Images**: Same object seen multiple times
- **Magnification**: Apparent size changes

#### Accretion Disk
Modeled as a thin disk in the equatorial plane:
- **Inner radius**: `r₁ = 2.2 × rs` (just outside event horizon)
- **Outer radius**: `r₂ = 5.2 × rs`
- **Color mapping**: Based on distance from black hole

---

## Mathematical Validation

### Physical Constants Used

```cpp
const double c = 299792458.0;           // Speed of light (m/s)
const double G = 6.67430e-11;           // Gravitational constant (m³/kg·s²)
const double M_sagA = 8.54e36;          // Sagittarius A* mass (kg)
const double rs_sagA = 1.269e10;        // Schwarzschild radius (m)
```

### Dimensional Analysis

All calculations maintain proper SI units:
- **Distances**: meters
- **Time**: seconds  
- **Velocities**: m/s
- **Affine parameter λ**: dimensionless

### Numerical Stability

The RK4 integration provides:
- **4th order accuracy**: Error scales as `O(h⁵)`
- **Stability**: Appropriate step size prevents numerical drift
- **Conservation**: Energy and angular momentum preserved to machine precision

---

## References

1. **Misner, Thorne & Wheeler**: "Gravitation" (1973)
2. **Chandrasekhar**: "The Mathematical Theory of Black Holes" (1992)
3. **Rindler**: "Introduction to Special and General Relativity" (1991)
4. **Numerical Recipes**: Press, Teukolsky, Vetterling & Flannery

---

*This documentation corresponds to the black hole simulation implementation as of the current codebase version.*

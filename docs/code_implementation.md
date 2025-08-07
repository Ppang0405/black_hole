# Code Implementation Guide

This document explains how the physics and mathematics are implemented in the actual code, with detailed examples and cross-references.

## Table of Contents

1. [Project Architecture](#project-architecture)
2. [Core Data Structures](#core-data-structures)
3. [Physics Implementation](#physics-implementation)
4. [Rendering Pipeline](#rendering-pipeline)
5. [GPU Compute Implementation](#gpu-compute-implementation)
6. [Camera and Controls](#camera-and-controls)
7. [Performance Optimizations](#performance-optimizations)
8. [Code Organization](#code-organization)

---

## Project Architecture

### File Structure and Dependencies

```
black_hole/
├── black_hole.cpp      # Main 3D simulation (CPU + GPU)
├── geodesic.comp       # GPU compute shader
├── 2D_lensing.cpp      # 2D educational version
├── CPU-geodesic.cpp    # CPU-only 3D version
├── ray_tracing.cpp     # Basic ray tracer demo
└── docs/               # Documentation
    ├── physics_and_mathematics.md
    ├── numerical_methods.md
    └── code_implementation.md
```

### Dependency Chain

```
OpenGL/GLFW/GLEW → Engine → Camera/Physics → Ray Tracing → Visualization
```

---

## Core Data Structures

### Ray Structure

The ray is the fundamental unit of computation, representing a light path through spacetime.

#### CPU Implementation (`CPU-geodesic.cpp`)

```cpp
struct Ray {
    // Cartesian coordinates
    double x, y, z;
    
    // Spherical coordinates  
    double r, theta, phi;
    
    // Velocity components in spherical basis
    double dr, dtheta, dphi;
    
    // Conserved quantities
    double E, L;  // Energy and angular momentum
    
    // Constructor: Initialize from position and direction
    Ray(vec3 pos, vec3 dir) : x(pos.x), y(pos.y), z(pos.z) {
        // Convert position to spherical coordinates
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z / r);
        phi = atan2(y, x);
        
        // Transform direction to spherical basis
        double dx = dir.x, dy = dir.y, dz = dir.z;
        dr     = sin(theta)*cos(phi)*dx + sin(theta)*sin(phi)*dy + cos(theta)*dz;
        dtheta = cos(theta)*cos(phi)*dx + cos(theta)*sin(phi)*dy - sin(theta)*dz;
        dtheta /= r;
        dphi   = -sin(phi)*dx + cos(phi)*dy;
        dphi  /= (r * sin(theta));
        
        // Calculate conserved quantities
        L = r * r * sin(theta) * dphi;
        double f = 1.0 - SagA.r_s / r;
        double dt_dλ = sqrt((dr*dr)/f + r*r*(dtheta*dtheta + sin(theta)*sin(theta)*dphi*dphi));
        E = f * dt_dλ;
    }
    
    void step(double dλ, double rs) {
        if (r <= rs) return;  // Stop at event horizon
        rk4Step(*this, dλ, rs);
        
        // Update Cartesian coordinates
        x = r * sin(theta) * cos(phi);
        y = r * sin(theta) * sin(phi);
        z = r * cos(theta);
    }
};
```

#### GPU Implementation (`geodesic.comp`)

```glsl
struct Ray {
    float x, y, z, r, theta, phi;
    float dr, dtheta, dphi;
    float E, L;
};

Ray initRay(vec3 pos, vec3 dir) {
    Ray ray;
    ray.x = pos.x; ray.y = pos.y; ray.z = pos.z;
    ray.r = length(pos);
    ray.theta = acos(pos.z / ray.r);
    ray.phi = atan(pos.y, pos.x);
    
    // Direction transformation (same logic as CPU)
    float dx = dir.x, dy = dir.y, dz = dir.z;
    ray.dr     = sin(ray.theta)*cos(ray.phi)*dx + sin(ray.theta)*sin(ray.phi)*dy + cos(ray.theta)*dz;
    ray.dtheta = (cos(ray.theta)*cos(ray.phi)*dx + cos(ray.theta)*sin(ray.phi)*dy - sin(ray.theta)*dz) / ray.r;
    ray.dphi   = (-sin(ray.phi)*dx + cos(ray.phi)*dy) / (ray.r * sin(ray.theta));
    
    // Conserved quantities
    ray.L = ray.r * ray.r * sin(ray.theta) * ray.dphi;
    float f = 1.0 - SagA_rs / ray.r;
    float dt_dL = sqrt((ray.dr*ray.dr)/f + ray.r*ray.r*(ray.dtheta*ray.dtheta + sin(ray.theta)*sin(ray.theta)*ray.dphi*ray.dphi));
    ray.E = f * dt_dL;
    
    return ray;
}
```

### Black Hole Structure

```cpp
struct BlackHole {
    vec3 position;    // Center position
    double mass;      // Mass in kg
    double radius;    // Physical radius (not used for event horizon)
    double r_s;       // Schwarzschild radius
    
    BlackHole(vec3 pos, float m) : position(pos), mass(m) {
        r_s = 2.0 * G * mass / (c*c);  // Event horizon radius
    }
    
    bool Intercept(float px, float py, float pz) const {
        double dx = double(px) - double(position.x);
        double dy = double(py) - double(position.y);
        double dz = double(pz) - double(position.z);
        double dist2 = dx * dx + dy * dy + dz * dz;
        return dist2 < r_s * r_s;  // Inside event horizon
    }
};

// Sagittarius A* black hole at galactic center
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36);
```

### Object Data Structure

For additional gravitating bodies and visual objects:

```cpp
struct ObjectData {
    vec4 posRadius;   // xyz = position, w = radius
    vec4 color;       // rgb = color, a = alpha
    float mass;       // Mass for gravitational effects
    vec3 velocity;    // Velocity for orbital mechanics
};

// Example objects in the simulation
vector<ObjectData> objects = {
    { vec4(4e11f, 0.0f, 0.0f, 4e10f), vec4(1,1,0,1), 1.98892e30 },  // Yellow star
    { vec4(0.0f, 0.0f, 4e11f, 4e10f), vec4(1,0,0,1), 1.98892e30 },  // Red star
    { vec4(0.0f, 0.0f, 0.0f, SagA.r_s), vec4(0,0,0,1), SagA.mass }   // Black hole
};
```

---

## Physics Implementation

### Geodesic Equations Implementation

The heart of the simulation is solving the geodesic equations in Schwarzschild spacetime.

#### Right-Hand Side Function

```cpp
void geodesicRHS(const Ray& ray, double rhs[6], double rs) {
    double r = ray.r;
    double theta = ray.theta;
    double dr = ray.dr;
    double dtheta = ray.dtheta;
    double dphi = ray.dphi;
    double E = ray.E;
    
    double f = 1.0 - rs / r;           // Lapse function
    double dt_dlambda = E / f;         // Time coordinate velocity
    
    // First derivatives (velocities)
    rhs[0] = dr;       // dr/dλ = dr
    rhs[1] = dtheta;   // dθ/dλ = dθ  
    rhs[2] = dphi;     // dφ/dλ = dφ
    
    // Second derivatives (accelerations) from Schwarzschild null geodesics
    
    // d²r/dλ² = gravitational + centrifugal terms
    rhs[3] = - (rs / (2 * r * r)) * f * dt_dlambda * dt_dlambda     // Gravitational attraction
           + (rs / (2 * r * r * f)) * dr * dr                      // Relativistic correction
           + r * (dtheta * dtheta + sin(theta) * sin(theta) * dphi * dphi);  // Centrifugal
    
    // d²θ/dλ² = Coriolis + centrifugal terms
    rhs[4] = - (2.0 / r) * dr * dtheta                             // Coriolis-type
           + sin(theta) * cos(theta) * dphi * dphi;                // Centrifugal (φ)
    
    // d²φ/dλ² = Coriolis terms
    rhs[5] = - (2.0 / r) * dr * dphi                               // Coriolis (r)
           - 2.0 * cos(theta) / sin(theta) * dtheta * dphi;        // Coriolis (θ)
}
```

#### RK4 Integration Step

```cpp
void rk4Step(Ray& ray, double dλ, double rs) {
    double y0[6] = { ray.r, ray.theta, ray.phi, ray.dr, ray.dtheta, ray.dphi };
    double k1[6], k2[6], k3[6], k4[6], temp[6];
    
    // k1: Slope at beginning of interval
    geodesicRHS(ray, k1, rs);
    
    // k2: Slope at midpoint using k1
    addState(y0, k1, dλ/2.0, temp);
    Ray r2 = ray; 
    r2.r = temp[0]; r2.theta = temp[1]; r2.phi = temp[2];
    r2.dr = temp[3]; r2.dtheta = temp[4]; r2.dphi = temp[5];
    geodesicRHS(r2, k2, rs);
    
    // k3: Slope at midpoint using k2
    addState(y0, k2, dλ/2.0, temp);
    Ray r3 = ray;
    r3.r = temp[0]; r3.theta = temp[1]; r3.phi = temp[2];
    r3.dr = temp[3]; r3.dtheta = temp[4]; r3.dphi = temp[5];
    geodesicRHS(r3, k3, rs);
    
    // k4: Slope at end using k3
    addState(y0, k3, dλ, temp);
    Ray r4 = ray;
    r4.r = temp[0]; r4.theta = temp[1]; r4.phi = temp[2];
    r4.dr = temp[3]; r4.dtheta = temp[4]; r4.dphi = temp[5];
    geodesicRHS(r4, k4, rs);
    
    // Weighted average of slopes
    ray.r      += (dλ/6.0)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    ray.theta  += (dλ/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    ray.phi    += (dλ/6.0)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
    ray.dr     += (dλ/6.0)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
    ray.dtheta += (dλ/6.0)*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4]);
    ray.dphi   += (dλ/6.0)*(k1[5] + 2*k2[5] + 2*k3[5] + k4[5]);
}
```

### Gravity Simulation

For the N-body gravitational dynamics of objects:

```cpp
// From black_hole.cpp main loop
for (auto& obj : objects) {
    for (auto& obj2 : objects) {
        if (&obj == &obj2) continue;  // Skip self-interaction
        
        float dx = obj2.posRadius.x - obj.posRadius.x;
        float dy = obj2.posRadius.y - obj.posRadius.y;
        float dz = obj2.posRadius.z - obj.posRadius.z;
        float distance = sqrt(dx * dx + dy * dy + dz * dz);
        
        if (distance > 0) {
            vector<double> direction = {dx/distance, dy/distance, dz/distance};
            double Gforce = (G * obj.mass * obj2.mass) / (distance * distance);
            double acc1 = Gforce / obj.mass;
            
            std::vector<double> acc = {direction[0]*acc1, direction[1]*acc1, direction[2]*acc1};
            
            if (Gravity) {  // Toggle gravity on/off
                obj.velocity.x += acc[0];
                obj.velocity.y += acc[1];
                obj.velocity.z += acc[2];
                
                obj.posRadius.x += obj.velocity.x;
                obj.posRadius.y += obj.velocity.y;
                obj.posRadius.z += obj.velocity.z;
            }
        }
    }
}
```

---

## Rendering Pipeline

### Engine Architecture

The `Engine` class manages the OpenGL rendering pipeline:

```cpp
struct Engine {
    GLFWwindow* window;
    GLuint quadVAO;           // Full-screen quad for texture rendering
    GLuint texture;           // Ray-traced result texture
    GLuint shaderProgram;     // Vertex/fragment shaders
    GLuint computeProgram;    // Compute shader for ray tracing
    
    // Uniform Buffer Objects for GPU data transfer
    GLuint cameraUBO;         // Camera parameters
    GLuint diskUBO;           // Accretion disk parameters
    GLuint objectsUBO;        // Object array data
    
    // Grid rendering for spacetime visualization
    GLuint gridVAO, gridVBO, gridEBO;
    int gridIndexCount;
    
    // Resolution settings
    int WIDTH = 800, HEIGHT = 600;           // Window resolution
    int COMPUTE_WIDTH = 200, COMPUTE_HEIGHT = 150;  // Compute resolution
};
```

### Ray Tracing Loop

The main simulation loop combines physics and rendering:

```cpp
int main() {
    setupCameraCallbacks(engine.window);
    
    while (!glfwWindowShouldClose(engine.window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Update object positions with gravity
        for (auto& obj : objects) {
            // ... gravity calculation (see above) ...
        }
        
        // Generate curved spacetime grid
        engine.generateGrid(objects);
        
        // Set up camera matrices
        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 proj = perspective(radians(60.0f), 
                               float(engine.COMPUTE_WIDTH)/engine.COMPUTE_HEIGHT, 
                               1e9f, 1e14f);
        mat4 viewProj = proj * view;
        
        // Render spacetime grid
        engine.drawGrid(viewProj);
        
        // Run GPU ray tracer
        engine.dispatchCompute(camera);
        
        // Display result
        engine.drawFullScreenQuad();
        
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
}
```

### Spacetime Grid Generation

Visualizing curved spacetime as a deformed grid:

```cpp
void generateGrid(const vector<ObjectData>& objects) {
    const int gridSize = 25;
    const float spacing = 1e10f;
    
    vector<vec3> vertices;
    vector<GLuint> indices;
    
    for (int z = 0; z <= gridSize; ++z) {
        for (int x = 0; x <= gridSize; ++x) {
            float worldX = (x - gridSize / 2) * spacing;
            float worldZ = (z - gridSize / 2) * spacing;
            float y = 0.0f;
            
            // Warp grid using Schwarzschild geometry
            for (const auto& obj : objects) {
                vec3 objPos = vec3(obj.posRadius);
                double mass = obj.mass;
                double r_s = 2.0 * G * mass / (c * c);
                double dx = worldX - objPos.x;
                double dz = worldZ - objPos.z;
                double dist = sqrt(dx * dx + dz * dz);
                
                if (dist > r_s) {
                    // Curvature formula: approximation of spacetime warping
                    double deltaY = 2.0 * sqrt(r_s * (dist - r_s));
                    y += static_cast<float>(deltaY) - 3e10f;
                } else {
                    // Deep pit inside event horizon
                    y += 2.0f * static_cast<float>(sqrt(r_s * r_s)) - 3e10f;
                }
            }
            
            vertices.emplace_back(worldX, y, worldZ);
        }
    }
    
    // Generate line indices for wireframe rendering
    for (int z = 0; z < gridSize; ++z) {
        for (int x = 0; x < gridSize; ++x) {
            int i = z * (gridSize + 1) + x;
            indices.push_back(i);
            indices.push_back(i + 1);           // Horizontal line
            indices.push_back(i);
            indices.push_back(i + gridSize + 1); // Vertical line
        }
    }
    
    // Upload to GPU
    glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vec3), vertices.data(), GL_DYNAMIC_DRAW);
    // ... bind indices and vertex attributes ...
}
```

---

## GPU Compute Implementation

### Compute Shader Main Function

The GPU ray tracer processes thousands of rays in parallel:

```glsl
void main() {
    int WIDTH  = cam.moving ? 200 : 200;   // Adaptive resolution
    int HEIGHT = cam.moving ? 150 : 150;
    
    ivec2 pix = ivec2(gl_GlobalInvocationID.xy);
    if (pix.x >= WIDTH || pix.y >= HEIGHT) return;
    
    // Generate ray direction for this pixel
    float u = (2.0 * (pix.x + 0.5) / WIDTH - 1.0) * cam.aspect * cam.tanHalfFov;
    float v = (1.0 - 2.0 * (pix.y + 0.5) / HEIGHT) * cam.tanHalfFov;
    vec3 dir = normalize(u * cam.camRight - v * cam.camUp + cam.camForward);
    
    // Initialize ray from camera
    Ray ray = initRay(cam.camPos, dir);
    
    vec4 color = vec4(0.0);
    vec3 prevPos = vec3(ray.x, ray.y, ray.z);
    
    bool hitBlackHole = false;
    bool hitDisk = false;
    bool hitObject = false;
    
    int steps = cam.moving ? 60000 : 60000;
    
    // Integrate ray through spacetime
    for (int i = 0; i < steps; ++i) {
        // Check termination conditions
        if (intercept(ray, SagA_rs)) { 
            hitBlackHole = true; 
            break; 
        }
        
        // Advance ray one step
        rk4Step(ray, D_LAMBDA);
        
        vec3 newPos = vec3(ray.x, ray.y, ray.z);
        
        // Check for disk crossing
        if (crossesEquatorialPlane(prevPos, newPos)) { 
            hitDisk = true; 
            break; 
        }
        
        // Check object intersections
        if (interceptObject(ray)) { 
            hitObject = true; 
            break; 
        }
        
        prevPos = newPos;
        
        // Escape condition
        if (ray.r > ESCAPE_R) break;
    }
    
    // Determine pixel color based on ray fate
    if (hitDisk) {
        double r = length(vec3(ray.x, ray.y, ray.z)) / disk_r2;
        vec3 diskColor = vec3(1.0, r, 0.2);  // Orange-red disk
        color = vec4(diskColor, r);
    } else if (hitBlackHole) {
        color = vec4(0.0, 0.0, 0.0, 1.0);    // Black hole
    } else if (hitObject) {
        // Simple shading for objects
        vec3 P = vec3(ray.x, ray.y, ray.z);
        vec3 N = normalize(P - hitCenter);
        vec3 V = normalize(cam.camPos - P);
        float ambient = 0.1;
        float diff = max(dot(N, V), 0.0);
        float intensity = ambient + (1.0 - ambient) * diff;
        vec3 shaded = objectColor.rgb * intensity;
        color = vec4(shaded, objectColor.a);
    } else {
        color = vec4(0.0);  // Background (space)
    }
    
    // Write result to output image
    imageStore(outImage, pix, color);
}
```

### Uniform Buffer Object Management

Efficient data transfer to GPU:

```cpp
void uploadCameraUBO(const Camera& cam) {
    struct UBOData {
        vec3 pos; float _pad0;
        vec3 right; float _pad1;  
        vec3 up; float _pad2;
        vec3 forward; float _pad3;
        float tanHalfFov;
        float aspect;
        bool moving;
        int _pad4;
    } data;
    
    // Calculate camera basis vectors
    vec3 fwd = normalize(cam.target - cam.position());
    vec3 up = vec3(0, 1, 0);
    vec3 right = normalize(cross(fwd, up));
    up = cross(right, fwd);
    
    // Fill data structure
    data.pos = cam.position();
    data.right = right;
    data.up = up;
    data.forward = fwd;
    data.tanHalfFov = tan(radians(60.0f * 0.5f));
    data.aspect = float(WIDTH) / float(HEIGHT);
    data.moving = cam.dragging || cam.panning;
    
    // Upload to GPU
    glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UBOData), &data);
}
```

---

## Camera and Controls

### Camera System Implementation

The camera orbits around the black hole center:

```cpp
struct Camera {
    vec3 target = vec3(0.0f, 0.0f, 0.0f);  // Always look at black hole
    float radius = 6.34194e10f;             // Distance from center
    float minRadius = 1e10f, maxRadius = 1e12f;
    
    float azimuth = 0.0f;                   // Horizontal rotation
    float elevation = M_PI / 2.0f;          // Vertical angle
    
    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = 25e9f;
    
    bool dragging = false;
    bool panning = false;
    bool moving = false;
    double lastX = 0.0, lastY = 0.0;
    
    // Calculate world space position
    vec3 position() const {
        float clampedElevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    
    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);
        
        if (dragging && !panning) {
            // Orbit around black hole
            azimuth += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        }
        
        lastX = x; lastY = y;
        update();
    }
    
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        radius = clamp(radius, minRadius, maxRadius);
        update();
    }
};
```

### Input Handling

```cpp
void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera);
    
    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseButton(button, action, mods, win);
    });
    
    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseMove(x, y);
    });
    
    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processScroll(xoffset, yoffset);
    });
    
    glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int scancode, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processKey(key, scancode, action, mods);
    });
}
```

---

## Performance Optimizations

### Adaptive Quality Rendering

```cpp
void dispatchCompute(const Camera& cam) {
    // Lower resolution during camera movement
    int cw = cam.moving ? COMPUTE_WIDTH : 200;
    int ch = cam.moving ? COMPUTE_HEIGHT : 150;
    
    // Reallocate texture if size changed
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, cw, ch, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    
    // Bind compute program and data
    glUseProgram(computeProgram);
    uploadCameraUBO(cam);
    uploadDiskUBO();
    uploadObjectsUBO(objects);
    
    // Bind output texture
    glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);
    
    // Dispatch compute grid
    GLuint groupsX = (GLuint)std::ceil(cw / 16.0f);
    GLuint groupsY = (GLuint)std::ceil(ch / 16.0f);
    glDispatchCompute(groupsX, groupsY, 1);
    
    // Wait for completion
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
}
```

### Memory Layout Optimization

```cpp
// Structure of Arrays layout for GPU efficiency
struct UBOData {
    int numObjects;
    float _pad0, _pad1, _pad2;        // std140 alignment padding
    vec4 posRadius[16];               // All positions together
    vec4 color[16];                   // All colors together  
    float mass[16];                   // All masses together
} data;
```

### Early Ray Termination

Multiple exit conditions prevent unnecessary computation:

```glsl
for (int i = 0; i < steps; ++i) {
    if (intercept(ray, SagA_rs)) { hitBlackHole = true; break; }  // Event horizon
    rk4Step(ray, D_LAMBDA);
    if (crossesEquatorialPlane(prevPos, newPos)) { hitDisk = true; break; }  // Disk
    if (interceptObject(ray)) { hitObject = true; break; }        // Objects
    if (ray.r > ESCAPE_R) break;                                  // Escape
}
```

---

## Code Organization

### Modular Design

Each file has a specific purpose:

#### `black_hole.cpp` - Main Application
- Engine initialization and management
- Camera system and input handling
- OpenGL setup and rendering pipeline
- Compute shader dispatch and UBO management
- Main simulation loop

#### `geodesic.comp` - GPU Ray Tracer
- Ray initialization and integration
- Intersection testing (black hole, disk, objects)
- Shading and color determination
- Parallel pixel processing

#### `2D_lensing.cpp` - Educational Demo
- Simplified 2D physics for learning
- Direct OpenGL rendering
- Interactive ray creation
- Visual trail rendering

#### `CPU-geodesic.cpp` - Reference Implementation
- Full CPU implementation for comparison
- Higher precision double arithmetic
- Detailed error checking
- Performance benchmarking

### Constants and Configuration

Physical and numerical constants are centralized:

```cpp
// Physics constants
double c = 299792458.0;           // Speed of light (m/s)
double G = 6.67430e-11;           // Gravitational constant
double SagA_mass = 8.54e36;       // Sagittarius A* mass (kg)

// Numerical parameters
const double D_LAMBDA = 1e7;      // Integration step size
const int MAX_STEPS = 60000;      // Maximum ray steps
const double ESCAPE_R = 1e30;     // Escape radius

// Rendering parameters
const int COMPUTE_WIDTH = 200;    // Compute texture width
const int COMPUTE_HEIGHT = 150;   // Compute texture height
```

### Error Handling

```cpp
// Shader compilation error checking
GLint success;
glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
if (!success) {
    GLint logLen;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
    std::vector<char> log(logLen);
    glGetShaderInfoLog(shader, logLen, nullptr, log.data());
    std::cerr << "Shader compile error:\n" << log.data() << "\n";
    exit(EXIT_FAILURE);
}

// Ray integration bounds checking
if (ray.r <= 0.5 * rs) {
    // Too close to singularity
    return;
}
```

---

## Building and Running

### Dependencies

```cpp
#include <GL/glew.h>          // OpenGL extension loading
#include <GLFW/glfw3.h>       // Window management
#include <glm/glm.hpp>        // Mathematics library
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
```

### Compilation

#### Linux (Ubuntu/Debian)
```bash
# Install dependencies
sudo apt-get update
sudo apt-get install libglfw3-dev libglew-dev libglm-dev build-essential

# Compile main simulation
g++ -std=c++11 black_hole.cpp -lglfw -lGLEW -lGL -o black_hole

# Compile 2D version
g++ -std=c++11 2D_lensing.cpp -lglfw -lGLEW -lGL -o 2D_lensing

# Compile CPU version
g++ -std=c++11 CPU-geodesic.cpp -lglfw -lGLEW -lGL -o CPU-geodesic
```

#### Linux (CentOS/RHEL/Fedora)
```bash
# CentOS/RHEL (with EPEL)
sudo yum install epel-release
sudo yum install glfw-devel glew-devel glm-devel gcc-c++

# Fedora
sudo dnf install glfw-devel glew-devel glm-devel gcc-c++

# Compilation (same as Ubuntu)
g++ -std=c++11 black_hole.cpp -lglfw -lGLEW -lGL -o black_hole
```

#### macOS
```bash
# Install Homebrew (if needed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install glfw glew glm

# Compile with framework linking
g++ -std=c++11 black_hole.cpp -lglfw -lGLEW -framework OpenGL -o black_hole

# Alternative with explicit paths (Apple Silicon)
g++ -std=c++11 black_hole.cpp \
    -I/opt/homebrew/include \
    -L/opt/homebrew/lib \
    -lglfw -lGLEW -framework OpenGL \
    -o black_hole
```

#### Windows (MSYS2)
```bash
# Install MSYS2 from https://www.msys2.org/
# In MSYS2 terminal:

# Update packages
pacman -Syu

# Install dependencies
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-glfw mingw-w64-x86_64-glew mingw-w64-x86_64-glm

# Compile with Windows libraries
g++ -std=c++11 black_hole.cpp -lglfw3 -lglew32 -lopengl32 -lgdi32 -o black_hole.exe
```

#### Cross-platform CMake
Create `CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.10)
project(BlackHoleSimulation)

set(CMAKE_CXX_STANDARD 11)

find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(black_hole black_hole.cpp)
target_link_libraries(black_hole glfw GLEW::GLEW OpenGL::GL)
```

Build:
```bash
mkdir build && cd build
cmake ..
make  # Linux/macOS
cmake --build .  # Windows
```

### Runtime Controls

- **Mouse Drag**: Orbit camera around black hole
- **Mouse Scroll**: Zoom in/out
- **G Key**: Toggle gravity simulation
- **Right Click**: Enable gravity (alternative)
- **ESC**: Exit application

---

This implementation guide provides a complete picture of how the mathematical physics is translated into working code, from low-level ray integration to high-level rendering and user interaction.

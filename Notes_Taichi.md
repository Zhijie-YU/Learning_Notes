# Taichi
[TOC]
## Taichi basics
The gene of Taichi is parallel computing.
### Data format
Tensor is a multidim array whose elements can be everything even matrices.
```python
import taichi as ti
ti.init()

a = ti.var(dt=ti.f32, shape=(42,63)) # a tensor of 42X63 scalars
b = ti.Vector(3, dt=ti.f32, shape=4) # a tensor of 4X3D vectors
c = ti.Matrix(2,2,dt=ti.f32,shape=(3,5)) # a tensor of 3X5 2X2 matrices
# ti.* defines the type of each element and the first part defines the size of each element
loss = ti.var(dt=ti.f32, shape=()) # this is a scalar defined in tensor form
loss[None]=3 # use this to assign value to a scalar in a tensor form
```

### Kernel and function
In Taichi, kernel is the computation function.
Kernels must be decorated with @ti.kernel. They can call functions but cannot be call other kernels.
Taichi functions can be called directly by Taichi kernels and other Taichi functions but not python. Only **one return** is supported up to now. They must be decorated with @ti.func.

Element-wise product *; marix product @.

### For loops
For loops in Taichi have 2 forms.
+ Range-for loops: Similar to Python. Will be parallelized when used at the outermost scope.
+ Struct-for loops: Iterates over (sparse) tensor elements.
For loops at the outermost scope in Taichi is **automatically parallelized**.

```py
ti.init(arch=ti.gpu)
n = 320
pixels = ti.var(dt=ti.f32, shape=(2*n, n))
@ti.kernel
def paint(t:ti.f32):
    for i,j in pixels:
        pixels[i,j] = i*3+j*4+t
```

### Atomic operations
Atomic operation is designed to deal with parallel computing. An atomic operation will go from start to finish without interruption of other threads.
In Taichi, augmented assignments (x[i]+=1)
```py
a[None] += 1 [right]
ti.atomic_add(a[None], 1) [right]
a[None] = a[None] + 1 [wrong]
```

### Scope
Taichi-scope: in @ti.kernel or @ti.func. Compiled in Taichi and run in parallel.

Python-scope: Compiled in Python.

### Phases of a Taichi program
+ Initialization: ti.init(...)
+ Tensor allocation: ti.var, ti.Vector, ti.Matrix
+ Computation (lauch kernels...)
+ Optional: restart the Taichi system (clear memory, destroy variables and kernels...) ti.reset()
  

> :grey_exclamation::eyes: After the first ti.kernel, no more tensor allocation is allowed.

```py
# fractal.py
import taichi as ti
ti.init(arch=ti.cpu)

n=320
pixels = ti.var(dt=ti.f32, shape=(2*n, n))

@ti.func
def complex_sqr(z):
    return ti.Vector([z[0]**2 - z[1]**2, z[1] * z[0] *2])

@ti.kernel
def paint(t: ti.f32):
    for i,j in pixels:
        c = ti.Vector([-0.8, ti.cos(t) * 0.2])
        z = ti.Vector([i/n - 1, j/n - 0.5]) * 2
        iterations = 0
        while z.norm() < 20 and iterations < 50:
            z = complex_sqr(z) + c
            iterations += 1
        pixels[i, j] = 1 - iterations * 0.02

gui = ti.GUI("Julia Set", res=(n * 2, n))
for i in range(1000000):
    paint(i * 0.03)
    gui.set_image(pixels)
    gui.show()
```

### Debug mode
debug = True (cpu only)
```py
ti.init(debug=True, arch=ti.cpu)
a = ti.var(dt=ti.f32, shape=(10))
b = ti.var(dt=ti.f32, shape=(10))
@ti.kernel
def shift():
    for i in range(10):
        a[i] = b[i+1]
# bound checker is only activated in debug mode to save time in normal mode.
shift()
```

## Lagrangian and Eulerian View

Lagrangian view: move with object.

Eulerian view: static.

### Lagrangian simulation approaches (1)

#### Mass-spring system
$$
\begin{aligned}
    \mathbf{f}{ij}&=-k(\|\mathbf{x}_i-\mathbf{x}_j\|_2-l_{ij})(\widehat{\mathbf{x}_i-\mathbf{x}_j})\quad(\rm Hooke's\,law)\\
    \mathbf{f}_i&=\sum_j^{j\neq i}\mathbf{f}_{ij}\\
    \frac{\partial \mathbf{v}_i}{\partial t}&=\frac{1}{m_i}\mathbf{f}_i\quad(\rm Newton's\,second\,law\,of\,motion)\\
    \frac{\partial \mathbf{x}_i}{\partial t}&=\mathbf{v}_i\\
\end{aligned}
$$

$\widehat{\mathbf{x}_i-\mathbf{x}_j}$: direction vector from particle $i$ to particle $j$ (unit vector).
$\widehat{\Box}$ means **normalization**.

**Implicit time integration:**
![](Taichi_images/ms_1.png)
with taylor expansion
![](Taichi_images/ms_2.png)
![](Taichi_images/ms_3.png)
To solver this linear system, there are many methods like Jacobi iteration/Gauss-Seidel iteration or conjugate gradients(共轭梯度), etc.
**Unifying explict and implicit:**
![](Taichi_images/ms_4.png)
**Solve faster**
For system with millions of mass points and springs,
+ Sparse matrices
+ Conjugate gradients
+ Preconditioning
+ Use position-based dynamics(PBD)
+ Also some faster approches like Fast mass-spring system solver("Fast simulation of mass-spring systems" in ACM Transactions)



#### Time integration
:one: Forward Euler (explicit)
$$
\begin{aligned}
\mathbf{v}_{t+1}&=\mathbf{v}_t+\Delta t\frac{\mathbf{f}_t}{m}\\
\mathbf{x}_{t+1}&=\mathbf{x}_t+\Delta t\mathbf{v}_t\\
\end{aligned}
$$

:two: Semi implicit Euler (aka. symplectic Euler, explicit)
$$
\begin{aligned}
\mathbf{v}_{t+1}&=\mathbf{v}_t+\Delta t\frac{\mathbf{f}_t}{m}\\
\mathbf{x}_{t+1}&=\mathbf{x}_t+\Delta t\mathbf{v}_{t+1}\\
\end{aligned}
$$

```py
# mass_spring.py

```

:three: Backward Euler (often with Newton's method, implicit)

#### Explicit v.s. implicit time integration
Explicit (forward Euler, symplectic Euler, RK, ...)
$$\Delta t \leq c\sqrt{\frac{m}{k}} \quad c\approx1$$ 

Implicit (backward Euler, middle-point, ...)

#### Lagrangian fluid simulation: Smoothed particle hydrodynamics(SPH)
**Courant-Friedrichs-Lewy(CFL) condition**
Another threshold
![](Taichi_images/ms_5.png)
**Accelerating SPH: Neighborhood search**
![](Taichi_images/ms_6.png)

#### Output mp4 and gif in taichi
ti.imwrite(img, filename)
ti video -f 24 or ti video -f 60
ti git -i input.mp4

Make sure ffmpeg installed!

### Lagrangian simulation approaches (2)

#### Basics of deformation, elasticity and FEM
**Deformation**
Deformation map $\phi$:$$\mathbf{x}_{\rm deformed}=\phi(\mathbf{x}_{\rm rest})$$

> This relates rest material position with deformed material position.

Deformation gradient $\mathbf{F}$:$$\mathbf{F}:=\frac{\partial\mathbf{x}_{\rm deformed}}{\partial\mathbf{x}_{\rm rest}}$$

> Deformation gradients are translational invariant.
> $\phi_1=\phi(\mathbf{x}_{\rm rest})$ and $\phi_2=\phi(\mathbf{x}_{\rm rest})+\rm c$ have the same $\mathbf{F}$.

Deform/rest volume ratio $J=\det(\mathbf{F})$

**Elasticity**
Hyperelasticity
whose stress-strain relationship is defined by **strain energy density function**.
$$\psi=\psi(\mathbf{F})$$

There are different measures of stress:
* THe First Piola-Kirchhoff stress tensor (PK1):$\mathbf{P}(\mathbf{F})=\frac{\partial\psi(\mathbf{F})}{\partial\mathbf{F}}$ (easy to compute but in rest space)
* Kirchhoff stress: $\boldsymbol{\tau}$
* Cauchy stress tensor: $\boldsymbol{\sigma}$ (symmetric)

Relationship: $\boldsymbol{\tau}=J\boldsymbol{\sigma}=\mathbf{P}\mathbf{F}^T\quad\mathbf{P}=J\boldsymbol{\sigma}\mathbf{F}^{-T}\quad \rm Traction\,\mathbf{t}=\boldsymbol{\sigma}^T\boldsymbol{n}$

* Young's modulus $E = \frac{\sigma}{\varepsilon}$
* Bulk modulus $K = -V\frac{dP}{dV}$
* Poisson's ratio $\nu\in[0,0.5)$
* Lame's first parameter $\mu$; Lame's second parameter $\lambda$ (aka. shear modulus $G$)

conversion formula:
$$K=\frac{E}{3(1-2\nu)}\qquad\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}\qquad\mu=\frac{E}{2(1+\nu)}$$

Popular hyperelastic material models
* Neo-Hookean
  * $\psi(\mathbf{F})=\frac{\mu}{2}\sum_i[(\mathbf{F}^T\mathbf{F})_{ii}-1]-\mu\log(J)+\frac{\lambda}{2}\log^2(J)$
  * $\mathbf{P}(\mathbf{F})=\frac{\partial\psi(\mathbf{F})}{\partial\mathbf{F}}=2\mu(\mathbf{F}-\mathbf{R})+\lambda(J-1)J\mathbf{F}^{-T}$
* (Fixed) Corotated

**FEM**
Linear tetrahedral FEM
The deformation map $\phi$ is affine and thus deformation gradient $\mathbf{F}$ is **constant** within a single tetrahedral element:
$$\rm \mathbf{x}_{deformed}=\mathbf{F}\mathbf{x}_{rest}+\mathbf{b}$$

For every element $e$, its elastic potential energy
$$U(e)=\int_e\psi(\mathbf{F}(\mathbf{x}))\mathbf{x}=V_e\psi(\mathbf{F}_e)$$

For explicit scheme (semi-implicit)
$$\mathbf{v}_{t+1,i}=\mathbf{v}_{t,i}+\Delta t\frac{\mathbf{f}_{t,i}}{m_i}$$

$$\mathbf{x}_{t+1,i}=\mathbf{x}_{t,i}+\Delta t\mathbf{v}_{t+1,i}$$

$$\mathbf{f}_{t,i}\equiv-\frac{\partial U}{\partial \mathbf{x}_i}=-\sum_e\frac{\partial U(e)}{\partial \mathbf{x}_i}=-\sum_e V_e\frac{\partial\psi(\mathbf{F}_e)}{\partial\mathbf{F}_e}\frac{\partial \mathbf{F}_e}{\partial \mathbf{x}_i}=-\sum_eV_e\mathbf{P}(\mathbf{F}_e)\frac{\partial\mathbf{F}_e}{\partial\mathbf{x}_i}$$

#### Taichi programming language advanced features
##### ODOP
Data-oriented programming (DOP)
Objective data-oriented programming (ODOP)
+ 3 important decorators
  + Use @ti.data_oriented to decorate class.
  + Use @ti.kernel to decorate class members functions that are Taichi kernels.
  + Use @ti.func to decorate class members functions that are Taichi functions.

##### Metaprogramming
+ Allow to pass almost anything to Taichi kernels
+ Improve run-time performance by moving run-time costs to compile time
+ Achieve dimensionality independence
+ Simplify the development of Taichi standard library
```py
@ti.kernel
def copy(x: ti.template(), y: ti.template(), c: ti.f32):
    for i in x:
        y[i] = x[i] + c
```
Variable aliasing
##### Differentiable programming
reverse-mode automatic differentiation (AutoDiff)
$$f(x)\Rightarrow\frac{\partial f(x)}{\partial x}$$
##### Visualization



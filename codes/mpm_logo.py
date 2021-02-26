import taichi as ti
import numpy as np
ti.init(arch=ti.cpu) # Try to run on GPU
quality = 2 # Use a larger value for higher-res simulations
n_particles, n_grid = 9000 * quality ** 2, 128 * quality
dx, inv_dx = 1 / n_grid, float(n_grid)
dt = 1e-4 / quality
p_vol, p_rho = (dx * 0.5)**2, 1
p_mass = p_vol * p_rho
E, nu = 0.1e4, 0.2 # Young's modulus and Poisson's ratio
mu_0, lambda_0 = E / (2 * (1 + nu)), E * nu / ((1+nu) * (1 - 2 * nu)) # Lame parameters
x = ti.Vector.field(2, dtype=float, shape=n_particles) # position
v = ti.Vector.field(2, dtype=float, shape=n_particles) # velocity
C = ti.Matrix.field(2, 2, dtype=float, shape=n_particles) # affine velocity field
F = ti.Matrix.field(2, 2, dtype=float, shape=n_particles) # deformation gradient
material = ti.field(dtype=int, shape=n_particles) # material id
Jp = ti.field(dtype=float, shape=n_particles) # plastic deformation
grid_v = ti.Vector.field(2, dtype=float, shape=(n_grid, n_grid)) # grid node momentum/velocity
grid_m = ti.field(dtype=float, shape=(n_grid, n_grid)) # grid node mass

grid_g = ti.Vector.field(2, dtype=float, shape=(n_grid, n_grid))

@ti.kernel
def substep():
  for i, j in grid_m:
    grid_v[i, j] = [0, 0]
    grid_m[i, j] = 0
  for p in x: # Particle state update and scatter to grid (P2G)
    base = (x[p] * inv_dx - 0.5).cast(int)
    fx = x[p] * inv_dx - base.cast(float)
    # Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
    w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
    F[p] = (ti.Matrix.identity(float, 2) + dt * C[p]) @ F[p] # deformation gradient update
    h = ti.exp(10 * (1.0 - Jp[p])) # Hardening coefficient: snow gets harder when compressed
    if material[p] == 1: # jelly, make it softer
      h = 0.3
    mu, la = mu_0 * h, lambda_0 * h
    if material[p] == 0 or material[p] == 1: # liquid
      mu = 0.0
    U, sig, V = ti.svd(F[p])
    J = 1.0
    for d in ti.static(range(2)):
      new_sig = sig[d, d]
      if material[p] == 2:  # Snow
        new_sig = min(max(sig[d, d], 1 - 2.5e-2), 1 + 4.5e-3)  # Plasticity
      Jp[p] *= sig[d, d] / new_sig
      sig[d, d] = new_sig
      J *= new_sig
    if material[p] == 0:  # Reset deformation gradient to avoid numerical instability
      F[p] = ti.Matrix.identity(float, 2) * ti.sqrt(J)
    elif material[p] == 2:
      F[p] = U @ sig @ V.transpose() # Reconstruct elastic deformation gradient after plasticity
    stress = 2 * mu * (F[p] - U @ V.transpose()) @ F[p].transpose() + ti.Matrix.identity(float, 2) * la * J * (J - 1)
    stress = (-dt * p_vol * 4 * inv_dx * inv_dx) * stress
    affine = stress + p_mass * C[p]
    for i, j in ti.static(ti.ndrange(3, 3)): # Loop over 3x3 grid node neighborhood
      offset = ti.Vector([i, j])
      dpos = (offset.cast(float) - fx) * dx
      weight = w[i][0] * w[j][1]
      grid_v[base + offset] += weight * (p_mass * v[p] + affine @ dpos)
      grid_m[base + offset] += weight * p_mass
  for i, j in grid_m:
    if grid_m[i, j] > 0: # No need for epsilon here
      grid_v[i, j] = (1 / grid_m[i, j]) * grid_v[i, j] # Momentum to velocity
      
      grid_v[i, j][1] -= dt * grid_g[i,j][1] # gravity
      
      if i < 3 and grid_v[i, j][0] < 0:          grid_v[i, j][0] = 0 # Boundary conditions
      if i > n_grid - 3 and grid_v[i, j][0] > 0: grid_v[i, j][0] = 0
      if j < 3 and grid_v[i, j][1] < 0:          grid_v[i, j][1] = 0
      if j > n_grid - 3 and grid_v[i, j][1] > 0: grid_v[i, j][1] = 0
  for p in x: # grid to particle (G2P)
    base = (x[p] * inv_dx - 0.5).cast(int)
    fx = x[p] * inv_dx - base.cast(float)
    w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1.0) ** 2, 0.5 * (fx - 0.5) ** 2]
    new_v = ti.Vector.zero(float, 2)
    new_C = ti.Matrix.zero(float, 2, 2)
    for i, j in ti.static(ti.ndrange(3, 3)): # loop over 3x3 grid node neighborhood
      dpos = ti.Vector([i, j]).cast(float) - fx
      g_v = grid_v[base + ti.Vector([i, j])]
      weight = w[i][0] * w[j][1]
      new_v += weight * g_v
      new_C += 4 * inv_dx * weight * g_v.outer_product(dpos)
    v[p], C[p] = new_v, new_C
    x[p] += dt * v[p] # advection

@ti.func
def radius(theta):
  r = 0.0
  t = 3.4
  b = 10.0
  a = 1/(ti.log(b)*pow(b,t*np.pi))
  th = t * np.pi
  r0 = (th - a * (pow(b,th)-1))/50
  if theta < th:
    r = (theta - a * (pow(b,theta)-1))/50
  else:
    r = r0
  return r
#group_size = n_particles // 3
@ti.kernel
def initialize():
  for i,j in grid_g:
    grid_g[i,j] = [0, (n_grid - j)/n_grid * 100]   
  
  center = [0.5,0.5]
  r1 = radius(1000)
  r0 = r1*1.2
  print(r1)
  for i in range(n_particles):
    #x[i] = [ti.random() * 0.2 + 0.3 + 0.10 * (i // group_size), ti.random() * 0.2 + 0.05 + 0.32 * (i // group_size)]
    #material[i] = i // group_size # 0: fluid 1: jelly 2: snow
    x[i] = [ti.random() * r0 * 2 + center[0] - r0, ti.random() * r0 *2 + center[1] - r0]
    r_i = x[i]-center
    r_in = r_i.norm()
    theta_i = ti.asin(r_i[1]/r_in)
    theta_1 = 0.0
    theta_2 = 0.0
    
    if r_in <= r1:
      if theta_i >= 0:
        if r_i[0] * r_i[1] >= 0:
          theta_1 = theta_i + np.pi
        else:
          theta_1 = 2*np.pi - theta_i
      else:
        if r_i[0] * r_i[1] >= 0:
          theta_1 = -theta_i
        else:
          theta_1 = np.pi + theta_i

      if theta_i >= 0:
        if r_i[0] * r_i[1] >= 0:
          theta_2 = theta_i
        else:
          theta_2 = np.pi - theta_i
      else:
        if r_i[0] * r_i[1] >= 0:
          theta_2 = np.pi - theta_i
        else:
          theta_2 = 2*np.pi + theta_i

      r11 = radius(theta_1)
      r22 = radius(theta_2)
      flag = 0
      if r11 > r22:
        flag = 1
      for j in range(5):
        if flag == 0:
          rr = radius(theta_1)
          if r_in <= rr:
            break
          else:
            flag = 1
            theta_1 = theta_1 + 2*np.pi
        else:
           rr = radius(theta_2)
           if r_in <= rr:
             break
           else:
             flag = 0
             theta_2 = theta_2 + 2*np.pi
      material[i] = flag

      v[i] = ti.Matrix([0, 0])
      F[i] = ti.Matrix([[1, 0], [0, 1]])
      Jp[i] = 1
initialize()
gui = ti.GUI("Taichi MLS-MPM-LOGO", res=512, background_color=0x000000)
while not gui.get_event(ti.GUI.ESCAPE, ti.GUI.EXIT):
  for s in range(int(2e-3 // dt)):
    substep()
  colors = np.array([0xFFFFFF, 0x6495ED], dtype=np.uint32)
  gui.circles(x.to_numpy(), radius=1.5, color=colors[material.to_numpy()])
  gui.show() # Change to gui.show(f'{frame:06d}.png') to write images to disk

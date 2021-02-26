import taichi as ti
from pyevtk.hl import pointsToVTK
import numpy as np
import time
import os

from DFSPH_SOLVER import SPHSolver

ti.init(arch=ti.cpu)

def main():
    dynamic_allocate = False
    adaptive_time_step = True

    sim_physical_time = 5.0
    max_frame = 50000

    aabb = (0, 0, 0, 5, 5, 8)
    dx = 0.1
    up, bottom, left, right, front, back = np.array([aabb[5], aabb[2], aabb[0], aabb[3], aabb[1], aabb[4]])
    orgX = aabb[0]
    orgY = aabb[1]
    orgZ = aabb[2]
    sizeX = aabb[3] - aabb[0]
    sizeY = aabb[4] - aabb[1]
    sizeZ = aabb[5] - aabb[2]
    thickness = 0.3

    sph = SPHSolver([up, bottom, left, right, front, back],
                    alpha=0.30,
                    dx=dx,
                    max_time=sim_physical_time,
                    max_steps=max_frame,
                    dynamic_allocate=dynamic_allocate,
                    adaptive_time_step=adaptive_time_step)

    # Add fluid particles
    sph.add_cube(lower_corner=[1.5, 1.5, 2],
                 cube_size=[2, 2, 4],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_fluid)
    
    # Add boundary
    sph.add_cube(lower_corner=[orgX-thickness, orgY-thickness, orgZ-thickness],
                 cube_size=[sizeX+2*thickness, sizeY+2*thickness, thickness],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_bound)
    
    sph.add_cube(lower_corner=[orgX-thickness, orgY-thickness, sizeZ],
                 cube_size=[sizeX+2*thickness, sizeY+2*thickness, thickness],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_bound)
    
    sph.add_cube(lower_corner=[orgX-thickness, orgY-thickness, orgZ],
                 cube_size=[sizeX+2*thickness, thickness, sizeZ],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_bound)
    
    sph.add_cube(lower_corner=[orgX-thickness, orgY+sizeY, orgZ],
                 cube_size=[sizeX+2*thickness, thickness, sizeZ],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_bound)
    
    sph.add_cube(lower_corner=[orgX-thickness, orgY, orgZ],
                 cube_size=[thickness, sizeY, sizeZ],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_bound)
    
    sph.add_cube(lower_corner=[orgX+sizeX, orgY, orgZ],
                 cube_size=[thickness, sizeY, sizeZ],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_bound)   

    t = 0.0
    frame = 0
    total_start = time.process_time()
    while frame < max_frame and t < sim_physical_time:
        dt = sph.step(frame, t, total_start)
        particles = sph.particle_info()
        
        if frame % 100 == 0:
            if 1:
                mat = particles["material"]
                fluid = [i for i,m in enumerate(mat) if m]
                
                pos = particles["position"][np.array(fluid)]
                vel = particles["velocity"][np.array(fluid)]
                
                density = particles["density"][np.array(fluid)]
                ddensity = particles["d_density"][np.array(fluid)]
                pressure = particles["pressure"][np.array(fluid)]
                pos_x = copyData(pos[:,0])[np.array(fluid)]
                pos_y = copyData(pos[:,1])[np.array(fluid)]
                pos_z = copyData(pos[:,2])[np.array(fluid)]
                vel_x = copyData(vel[:,0])[np.array(fluid)]
                vel_y = copyData(vel[:,1])[np.array(fluid)]
                vel_z = copyData(vel[:,2])[np.array(fluid)]
                vel = np.linalg.norm(vel, axis=1)[np.array(fluid)]               
    
                pointsToVTK(f'./vtkData/frame_{frame:06d}', pos_x, pos_y, pos_z, 
                            data={"vel_x": vel_x, "vel_y": vel_y, "vel_z": vel_z, "vel": vel, "material": mat,
                                  "density": density, "ddensity": ddensity, "pressure": pressure})
            
        frame += 1
        t += dt

    print('Finish')

def copyData(value):
    num = value.size
    rvalue = np.zeros(num)
    for i in range(num):
        rvalue[i] = value[i]
    return rvalue
        
        

if __name__ == '__main__':
    main()
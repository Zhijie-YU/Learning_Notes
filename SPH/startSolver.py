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

    aabb = (-1, -1, 0, 10, 3, 8)
    dx = 0.1
    up, bottom, left, right, front, back = np.array([aabb[5], aabb[2], aabb[0], aabb[3], aabb[1], aabb[4]])

    sph = SPHSolver([up, bottom, left, right, front, back],
                    alpha=0.30,
                    dx=dx,
                    max_time=sim_physical_time,
                    max_steps=max_frame,
                    dynamic_allocate=dynamic_allocate,
                    adaptive_time_step=adaptive_time_step)

    # Add fluid particles
    sph.add_cube(lower_corner=[0, 0, 0],
                 cube_size=[2, 2, 4],
                 velocity=[0.0, 0.0, 0.0],
                 density=[1000],
                 material=SPHSolver.material_fluid)
    
    #sph.add_cube(lower_corner=[8, 0, 0],
    #             cube_size=[1, 1, 6],
    #             velocity=[0.0, 0.0, 0.0],
    #             density=[1000],
    #             material=SPHSolver.material_fluid)

    t = 0.0
    frame = 0
    total_start = time.process_time()
    while frame < max_frame and t < sim_physical_time:
        dt = sph.step(frame, t, total_start)
        particles = sph.particle_info()
        
        if frame % 100 == 0:
            if 1:
                pos = particles["position"]
                vel = particles["velocity"]
                mat = particles["material"]
                density = particles["density"]
                ddensity = particles["d_density"]
                pressure = particles["pressure"]
                pos_x = copyData(pos[:,0])
                pos_y = copyData(pos[:,1])
                pos_z = copyData(pos[:,2])
                vel_x = copyData(vel[:,0])
                vel_y = copyData(vel[:,1])
                vel_z = copyData(vel[:,2])
                vel = np.linalg.norm(vel, axis=1)
    
                pointsToVTK(f'./vtkData/frame_{frame:06d}', pos_x, pos_y, pos_z, 
                            data={"vel_x": vel_x, "vel_y": vel_y, "vel_z": vel_z, "vel": vel, "material": mat,
                                  "density": density, "ddensity": ddensity, "pressure": pressure})
            
        frame += 1
        t += dt
        #save_cnt += dt

    print('done')

def copyData(value):
    num = value.size
    rvalue = np.zeros(num)
    for i in range(num):
        rvalue[i] = value[i]
    return rvalue
        
        

if __name__ == '__main__':
    main()
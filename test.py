'''
# normal_fractal.py
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
'''
'''
# variant_fractal.py
import taichi as ti
ti.init(arch=ti.gpu)
n = 640
pixels = ti.var(dt=ti.f32, shape=(n * 2, n))

@ti.func
def complex_power(z, power: ti.i32):
    r = ti.sqrt(z[0]**2 + z[1]**2)
    theta = ti.atan2(z[1], z[0])
    return ti.Vector([r**power * ti.cos(power*theta), r**power * ti.sin(power*theta)])

@ti.kernel
def paint(t: ti.f32, power: ti.i32):
    for i, j in pixels:  # Parallized over all pixels
        # Julia Set
        freq = 1.0 / power
        c = ti.Vector([0.7885 * ti.cos(freq*t), 0.7885 * ti.sin(freq*t)])
        z = ti.Vector([i / n - 1, j / n - 0.5]) * 2

        iterations = 0
        while z.norm() < 20 and iterations < 50:
            z = complex_power(z, power) + c
            iterations += 1
        pixels[i, j] = 1 - iterations * 0.02

power = eval(input("Power of z -> "))
gui = ti.GUI("Julia Set", res=(n * 2, n))

for i in range(1000):
    paint(i * 0.03, power)
    gui.set_image(pixels)
    gui.show()


# mass_spring.py

import taichi as ti
tensor = ti.var(ti.f32 , shape =(4, 8, 16, 32, 64))
@ti.kernel
def print_tensor_size(x: ti.template ()):
    print(x.dim())
    for i in ti.static(range(x.dim())):
        print(x.shape ()[i])
print_tensor_size(tensor)
'''
'''
x1 = -2
y1 = 0
x2 = 3
y2 = -4
x3 = 0.05
y3 = 0.05
y0 = -3
# 求两点斜率
def grad(x1,y1,x2,y2):
    if x1 == x2:
        # 连线垂直于y=y0
        return 1e12
    else:
        return (y1-y2)/(x1-x2)

s1 = (y3 - y0) * (y1 - y0)
s2 = (y3 - y0) * (y2 - y0)
if s1 * s2 < 0:
    # 两个点位于y=y0两侧
    if s1 < 0:
        g = grad(x3,y3,x1,y1)
    else:
        g = grad(x3,y3,x2,y2)
else:
    # 两个点位于y=y0同侧
    g1 = grad(x3,y3,x1,y1)
    g2 = grad(x3,y3,x2,y2)
    if abs(g1) > abs(g2):
        g = g1
    else:
        g = g2
x0 = (y0 - y3)/g + x3
print(x0)
'''
'''
import numpy as np
a = np.zeros((3,3))
for i in range(3):
    if a[i,2] == 0:
        a[i,2] = 1
        print(type(a[0,0]))
print(a)

test = np.array([[-3, 1, 2], [4, 0, 2], [4, 2, 0]], dtype=np.float64)

for i in range(3):
    if test[i, 2] == 0:
        test[i,2] = 0.1
print(test)

print(np.tan(75*np.pi/180))
'''
'''
import taichi as ti
import numpy as np

ti.init(arch=ti.cpu)
dim = 2
n_particles = 2
x = ti.Vector(dim, dt=ti.f32, shape=n_particles)

@ti.kernel
def reset(mode: ti.i32):
    for i in range(n_particles):
        x[i] = [ti.random() * 0.6 + 0.2, ti.random() * 0.6 + 0.2]
    print(ti.random())


reset(1)
'''

import taichi as ti

ti.init()

pixels = ti.var(ti.u8, shape=(512, 512, 3))

@ti.kernel
def paint():
    for i, j, k in pixels:
        pixels[i, j, k] = ti.random() * 255

result_dir = "./results"
video_manager = ti.VideoManager(output_dir=result_dir, framerate=24, automatic_build=False)

for i in range(50):
    paint()

    pixels_img = pixels.to_numpy()
    video_manager.write_frame(pixels_img)
    print(f'\rFrame {i+1}/50 is recorded', end='')

print()
print('Exporting .mp4 and .gif videos...')
video_manager.make_video(gif=True, mp4=True)
print(f'MP4 video is saved to {video_manager.get_output_filename(".mp4")}')
print(f'GIF video is saved to {video_manager.get_output_filename(".gif")}')
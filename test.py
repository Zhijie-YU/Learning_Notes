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
'''

# mass_spring.py

import taichi as ti
tensor = ti.var(ti.f32 , shape =(4, 8, 16, 32, 64))
@ti.kernel
def print_tensor_size(x: ti.template ()):
    print(x.dim())
    for i in ti.static(range(x.dim())):
        print(x.shape ()[i])
print_tensor_size(tensor)


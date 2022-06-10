import numpy as np
from manim import *

class Questao1(ThreeDScene):
    def surface_func(self, u, v):
        return np.array([(2*np.cos(u) - 1)*np.cos(u),
                         (2*np.cos(u) - 1)*np.sin(u),
                         v])
    def plane_func(self, u, v, x, y):
        return self.surface_func(u, v) \
               + x * np.array([np.sin(u)*(1 - 4*np.cos(u)),
                               2*np.cos(2*u) - np.cos(u),
                               0]) \
               + y * np.array([0,
                               0,
                               1])

    def construct(self):
        axes = ThreeDAxes(x_range=[-7, 7], x_length=8)
        label_z = axes.get_z_axis_label(Tex("$z$"))
        label_y = axes.get_y_axis_label(Tex("$y$"))
        label_x = axes.get_x_axis_label(Tex("$x$"))
        labels = (label_z, label_y, label_x)

        surface = Surface(
            lambda u, v: axes.c2p(*self.surface_func(u, v)),
            u_range=[-PI, PI],
            v_range=[-0, 1],
            resolution=100
        )
        self.set_camera_orientation(theta=20*DEGREES, phi=45*DEGREES, zoom=1)

        self.add(surface, axes, *labels)

        u, v = PI/3, 0.5
        plane = Surface(
            lambda x, y: axes.c2p(*self.plane_func(u, v, x, y)),
            u_range=[-1, 1],
            v_range=[-1, 1],
            checkerboard_colors=['#ABCA29', '#6B8E23'],
            fill_opacity=.7
        )
        self.add(plane)


class Questao3(ThreeDScene):
    def surface_func(self, u, v):
        return np.array([v*np.cos(u),
                         v*np.sin(u),
                         u])
    def transform(self, x, y, z):
        if bool(np.isclose(np.cos(z), 0, atol=0.1)):
            return np.array([0, y / np.sin(z), z])
        else:
            return np.array([0, x / np.cos(z), z])

    def construct(self):
        axes = ThreeDAxes(x_range=[-7, 7], x_length=8)
        label_z = axes.get_z_axis_label(Tex("$z$"))
        label_y = axes.get_y_axis_label(Tex("$y$"))
        label_x = axes.get_x_axis_label(Tex("$x$"))
        labels = (label_z, label_y, label_x)

        surface = Surface(
            lambda u, v: axes.c2p(*self.surface_func(u, v)),
            u_range=[-1.5, 1.5],
            v_range=[-4, 4],
            resolution=75
        )

        self.set_camera_orientation(theta=45*DEGREES, phi=45*DEGREES)
        self.add(surface, axes, *labels)
        self.wait(1)
        self.play(
            ApplyPointwiseFunction(
                lambda p: self.transform(p[0], p[1], p[2]),
                surface,
                run_time=4
            )
        )
        self.wait(1)

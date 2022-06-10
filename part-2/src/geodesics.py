import numpy as np
from manim import *

class MyCylinder(ThreeDScene):
    def cylinder_func(self, u, v):
        return np.array([
            np.cos(u),
            np.sin(u),
            v,
        ])
    def geodesic_func(self, p_0, d, t):
        p = p_0 + t * d
        return self.cylinder_func(p[0], p[1])

    def construct(self):
        axes = ThreeDAxes(x_range=[-7, 7], x_length=8)
        label_z = axes.get_z_axis_label(Tex("$z$"))
        label_y = axes.get_y_axis_label(Tex("$y$"))
        label_x = axes.get_x_axis_label(Tex("$x$"))
        labels = (label_z, label_y, label_x)

        # cilinder = Surface(
        #     lambda u, v: axes.c2p(*self.cilinder_func(u, v)),
        #     u_range=[-PI, PI],
        #     v_range=[-5, 5],
        #     resolution=50
        # )
        cylinder = Cylinder(radius=1, height=5)
        self.set_camera_orientation(theta=270*DEGREES, phi=75*DEGREES)

        self.add(cylinder, axes, *labels)

        p_0 = np.array([0.0, 0.0])
        for d in np.array([[np.sqrt(3)/2, 1/2], [0.0, 0.7], [1.0, 0.0]]):
            geodesic = ParametricFunction(lambda t: \
                                          self.geodesic_func(p_0, d, t),
                                          t_range = np.array([-3, 3]),
                                          fill_opacity=0,
                                          stroke_opacity=.7,
                                          dt=1e-10).set_color(RED)
            self.add(geodesic)
        # self.begin_ambient_camera_rotation(rate=.6, about="theta")
        # self.wait(8)

class MySphere(ThreeDScene):
    def construct(self):
        def great_arc(t):
            return np.array([np.cos(t), np.sin(t), 0])

        def R_x(theta):
            return np.array([[1, 0, 0],
                             [0, np.cos(theta), -np.sin(theta)],
                             [0, np.sin(theta), np.cos(theta)]])
        def R_y(theta):
            return np.array([[np.cos(theta), 0, np.sin(theta)],
                             [0, 1, 0],
                             [-np.sin(theta), 0, np.cos(theta)]])
        def R_z(theta):
            return np.array([[np.cos(theta), -np.sin(theta), 0],
                             [np.sin(theta), np.cos(theta), 0],
                             [0, 0, 1]])

        axes = ThreeDAxes(x_range=[-7, 7], x_length=8)
        label_z = axes.get_z_axis_label(Tex("$z$"))
        label_y = axes.get_y_axis_label(Tex("$y$"))
        label_x = axes.get_x_axis_label(Tex("$x$"))
        labels = (label_z, label_y, label_x)

        sphere = Sphere(radius=1)
        self.set_camera_orientation(theta=90*DEGREES, phi=60*DEGREES)

        self.add(sphere, axes, *labels)

        for theta in np.array([PI/4, 3*PI/4, 5*PI/4, 7*PI/4]):
            R = R_x(theta) @ R_y(theta) @ R_z(theta)
            geodesic = ParametricFunction(lambda t: \
                                          R @ great_arc(t),
                                          t_range = np.array([0, TAU]),
                                          fill_opacity=0,
                                          stroke_opacity=.6,
                                          dt=1e-10).set_color(RED)
            self.add(geodesic)
        # self.begin_ambient_camera_rotation(rate=.6, about="theta")
        # self.wait(8)

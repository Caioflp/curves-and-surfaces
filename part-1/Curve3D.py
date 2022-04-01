import numpy as np
from manim import *
import sympy as sp

class Curve3D():

    def __init__(self, name: str, t: sp.Symbol, vector_sym, init_t: float, color: str,
                 speed=1, scale=1, offset=[0, 0, 0]):

        ## Symbolic attributes

        # Parameter
        self.t_sym = t * speed

        # Parametrized Vector
        self.x_sym = vector_sym[0].subs(t, self.t_sym) * scale + offset[0]
        self.y_sym = vector_sym[1].subs(t, self.t_sym) * scale + offset[1]
        self.z_sym = vector_sym[2].subs(t, self.t_sym) * scale + offset[2]

        # Velocity
        self.dx_sym = sp.diff(self.x_sym)
        self.dy_sym = sp.diff(self.y_sym)
        self.dz_sym = sp.diff(self.z_sym)

        # Acceleration
        self.ddx_sym = sp.diff(self.dx_sym)
        self.ddy_sym = sp.diff(self.dy_sym)
        self.ddz_sym = sp.diff(self.dz_sym)

        ## Numeric functions (used internaly)

        self.x_ = sp.lambdify(t, self.x_sym, "numpy")
        self.y_ = sp.lambdify(t, self.y_sym, "numpy")
        self.z_ = sp.lambdify(t, self.z_sym, "numpy")

        self.dx_ = sp.lambdify(t, self.dx_sym, "numpy")
        self.dy_ = sp.lambdify(t, self.dy_sym, "numpy")
        self.dz_ = sp.lambdify(t, self.dz_sym, "numpy")

        self.ddx_ = sp.lambdify(t, self.ddx_sym, "numpy")
        self.ddy_ = sp.lambdify(t, self.ddy_sym, "numpy")
        self.ddz_ = sp.lambdify(t, self.ddz_sym, "numpy")

        ## Manim Objects

        # Time - used to trace the curve, equals actual time_manim passed since
        #        the start of the animation
        self.time_m = ValueTracker(init_t)

        ## TEMPORARY SOLUTION
        ## Vector updaters SHOULD NOT HAVE TO BE HERE
        def time_m_updater(time_m, dt):
            time_m.increment_value(dt)
            self.velocity_m.update()
            self.normal_m.update()
            self.binormal_m.update()

        self.time_m.add_updater(time_m_updater)


        # Dot - mark the current location of the curve with a dot
        self.dot_m = Dot3D(point=self.vector(init_t), color=color)

        def dot_m_updater(dot_m):
            t = self.time_m.get_value()
            dot_m.move_to(self.vector(t))

        self.dot_m.add_updater(dot_m_updater)


        # Velocity vector on screen
        self.velocity_m = Arrow3D(start=self.vector(init_t),
                                  end=self.vector(init_t)+self.velocity(init_t),
                                  color=color)

        def velocity_m_updater(velocity_m):
            t = self.time_m.get_value()
            start = self.dot_m.get_center()
            end = start + self.velocity(t)
            velocity_m.set_start_and_end_attrs(start, end)

        self.velocity_m.add_updater(velocity_m_updater)


        # Acceleration vector on screen
        self.acceleration_m = Arrow3D(start=self.vector(init_t),
                                  end=self.vector(init_t)+self.acceleration(init_t),
                                  color=color)

        def acceleration_m_updater(acceleration_m):
            t = self.time_m.get_value()
            start = self.dot_m.get_center()
            end = start + self.acceleration(t)
            acceleration_m.set_start_and_end_attrs(start, end)

        self.acceleration_m.add_updater(acceleration_m_updater)


        # Normal (normalized acceleration) vector on screen
        self.normal_m = Arrow3D(start=self.vector(init_t),
                                  end=self.vector(init_t)+self.normal(init_t),
                                  color=color)

        def normal_m_updater(normal_m):
            t = self.time_m.get_value()
            start = self.dot_m.get_center()
            end = start + self.normal(t)
            normal_m.set_start_and_end_attrs(start, end)

        self.normal_m.add_updater(normal_m_updater)


        # Binormal
        self.binormal_m = Arrow3D(start=self.vector(init_t),
                                  end=self.vector(init_t)+self.binormal(init_t),
                                  color=color)

        def binormal_m_updater(binormal_m):
            t = self.time_m.get_value()
            start = self.dot_m.get_center()
            end = start + self.binormal(t)
            binormal_m.set_start_and_end_attrs(start, end)

        self.binormal_m.add_updater(binormal_m_updater)


        # Velocity label
        self.velocity_name_m = MathTex("\\" + name + r"^{\prime}").scale(.75)

        def velocity_name_m_updater(velocity_name_m):
            t = self.time_m.get_value()
            velocity_name_m.move_to(self.vector(t) + self.velocity(t)*1.2)

        self.velocity_name_m.add_updater(velocity_name_m_updater)


        # Acceleration label
        self.acceleration_name_m = MathTex("\\" + name + r"^{\prime\prime}").scale(.75)

        def acceleration_name_m_updater(acceleration_name_m):
            t = self.time_m.get_value()
            acceleration_name_m.move_to(self.vector(t) + self.acceleration(t)*1.2)

        self.acceleration_name_m.add_updater(acceleration_name_m_updater)


        # Normal label
        self.normal_name_m = MathTex("N_" + "\\" + name).scale(.75)

        def normal_name_m_updater(normal_name_m):
            t = self.time_m.get_value()
            normal_name_m.move_to(self.vector(t) + self.normal(t)*1.2)

        self.normal_name_m.add_updater(normal_name_m_updater)


        # Binormal label
        self.binormal_name_m = MathTex("B_" + "\\" + name).scale(.75)

        def binormal_name_m_updater(binormal_name_m):
            t = self.time_m.get_value()
            binormal_name_m.move_to(self.vector(t) + self.binormal(t)*1.2)

        self.binormal_name_m.add_updater(binormal_name_m_updater)


        # Path to be traced
        self.path_m = TracedPath(self.dot_m.get_center, stroke_width=3,
                                 stroke_color=color)


    def x(self, t):
        return self.x_(t)

    def y(self, t):
        return self.y_(t)

    def z(self, t):
        return self.z_(t)

    def vector(self, t):
        return np.array([self.x(t), self.y(t), self.z(t)])

    def dx(self, t):
        return self.dx_(t)

    def dy(self, t):
        return self.dy_(t)

    def dz(self, t):
        return self.dz_(t)

    def velocity(self, t):
        return np.array([self.dx(t), self.dy(t), self.dz(t)])

    def ddx(self, t):
        return self.ddx_(t)

    def ddy(self, t):
        return self.ddy_(t)

    def ddz(self, t):
        return self.ddz_(t)

    def acceleration(self, t):
        return np.array([self.ddx(t), self.ddy(t), self.ddz(t)])

    def normal(self, t):
        T_prime = self.acceleration(t)
        return T_prime / np.linalg.norm(T_prime)

    def binormal(self, t):
        return np.cross(self.velocity(t), self.normal(t))

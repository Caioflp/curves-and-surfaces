import numpy as np
from manim import *
import sympy as sp

class Curve2D():

    def __init__(self, name: str, t: sp.Symbol, vector_sym, init_t: float, color: str,
                 speed=1, scale=1, offset=[0, 0]):

        ## Symbolic attributes

        # Parameter
        self.t_sym = t * speed

        # Parametrized Vector
        self.x_sym = vector_sym[0].subs(t, self.t_sym) * scale + offset[0]
        self.y_sym = vector_sym[1].subs(t, self.t_sym) * scale + offset[1]

        # Velocity
        self.dx_sym = sp.diff(self.x_sym)
        self.dy_sym = sp.diff(self.y_sym)

        # Acceleration
        self.ddx_sym = sp.diff(self.dx_sym)
        self.ddy_sym = sp.diff(self.dy_sym)

        ## Numeric functions (used internaly)

        self.x_ = sp.lambdify(t, self.x_sym, "numpy")
        self.y_ = sp.lambdify(t, self.y_sym, "numpy")
        self.dx_ = sp.lambdify(t, self.dx_sym, "numpy")
        self.dy_ = sp.lambdify(t, self.dy_sym, "numpy")
        self.ddx_ = sp.lambdify(t, self.ddx_sym, "numpy")
        self.ddy_ = sp.lambdify(t, self.ddy_sym, "numpy")

        ## Manim Objects

        # Time - used to trace the curve, equals actual time_manim passed since
        #        the start of the animation
        self.time_m = ValueTracker(init_t)

        def time_m_updater(time_m, dt):
            time_m.increment_value(dt)

        self.time_m.add_updater(time_m_updater)


        # Dot - mark the current location of the curve with a dot
        self.dot_m = Dot(point=self.vector(init_t), color=color)

        def dot_m_updater(dot_m):
            t = self.time_m.get_value()
            dot_m.move_to(self.vector(t))

        self.dot_m.add_updater(dot_m_updater)


        # Velocity vector on screen
        self.velocity_m = Vector(self.velocity(init_t), color=color)

        def velocity_m_updater(velocity_m):
            t = self.time_m.get_value()
            start = self.dot_m.get_center()
            end = start + self.velocity(t)
            velocity_m.put_start_and_end_on(start, end)

        self.velocity_m.add_updater(velocity_m_updater)


        # Acceleration vector on screen
        self.acceleration_m = Vector(self.acceleration(init_t), color=color)

        def acceleration_m_updater(acceleration_m):
            t = self.time_m.get_value()
            start = self.dot_m.get_center()
            end = start + self.acceleration(t)
            acceleration_m.put_start_and_end_on(start, end)

        self.acceleration_m.add_updater(acceleration_m_updater)

        # Path to be traced
        self.path_m = TracedPath(self.dot_m.get_center, stroke_width=3,
                                 stroke_color=color)

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


    def x(self, t):
        return self.x_(t)

    def y(self, t):
        return self.y_(t)

    def vector(self, t):
        return np.array([self.x(t), self.y(t), 0])

    def dx(self, t):
        return self.dx_(t)

    def dy(self, t):
        return self.dy_(t)

    def velocity(self, t):
        return np.array([self.dx(t), self.dy(t), 0])

    def ddx(self, t):
        return self.ddx_(t)

    def ddy(self, t):
        return self.ddy_(t)

    def acceleration(self, t):
        return np.array([self.ddx(t), self.ddy(t), 0])

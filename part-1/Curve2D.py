import numpy as np
from manim import *
from sympy import *

class Curve2D():

    def __init__(self, name: str, t: Symbol, vector_sym, init_t: float, color: str,
                 speed=1, scale=1, offset=[0, 0]):

        ## Symbolic attributes

        # Parameter
        self.t_sym = t * speed

        # Parametrized Vector
        self.x_sym = vector_sym[0] * scale + offset[0]
        self.y_sym = vector_sym[1] * scale + offset[1]

        # Velocity
        self.dx_sym = diff(self.x_sym)
        self.dy_sym = diff(self.y_sym)

        # Acceleration
        self.ddx_sym = diff(self.dx_sym)
        self.ddy_sym = diff(self.dy_sym)

        ## Numeric functions (used internaly)

        self.x_ = lambdify(t, self.x_sym, "numpy")
        self.y_ = lambdify(t, self.y_sym, "numpy")
        self.dx_ = lambdify(t, self.dx_sym, "numpy")
        self.dy_ = lambdify(t, self.dy_sym, "numpy")
        self.ddx_ = lambdify(t, self.ddx_sym, "numpy")
        self.ddy_ = lambdify(t, self.ddy_sym, "numpy")

        ## Manim Objects

        # Time - used to trace the curve, equals actual time passed since
        #        the start of the animation
        self.time = ValueTracker(init_t)
        self.time.add_updater(lambda mobject, dt: mobject.increment_value(dt))

        # Dot - mark the current location of the curve with a dot
        self.dot_manim = Dot(point=self.vector(init_t), color=color)
        self.dot_manim.add_updater(lambda x: x.move_to(self.vector(self.time.get_value())))

        # Velocity vector on screen
        self.velocity_manim = Vector(self.velocity(init_t), color=color)
        self.velocity_manim.add_updater(lambda x: x.put_start_and_end_on(self.dot_manim.get_center(), self.dot_manim.get_center() + self.velocity(self.time.get_value())))

        # Acceleration vector on screen
        self.acceleration_manim = Vector(self.acceleration(init_t), color=color)
        self.acceleration_manim.add_updater(lambda x: x.put_start_and_end_on(self.dot_manim.get_center(), self.dot_manim.get_center() + self.acceleration(self.time.get_value())))

        # Path to be traced
        self.path_manim = TracedPath(self.dot_manim.get_center, stroke_width=3,
                                     stroke_color=color)

        # Velocity and acceleration
        self.velocity_name_manim = MathTex("\\" + name + r"^{\prime}").scale(.75)
        self.acceleration_name_manim = MathTex("\\" + name + r"^{\prime \prime}").scale(.75)

        self.velocity_name_manim.add_updater(lambda x: x.move_to(self.vector(self.time.get_value()) + self.velocity(self.time.get_value())*1.2))
        self.acceleration_name_manim.add_updater(lambda x: x.move_to(self.vector(self.time.get_value()) + self.acceleration(self.time.get_value())*1.2))

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

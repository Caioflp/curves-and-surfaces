import numpy as np
import manim
from manim import *
from sympy import *
from Curve2D import *

## Fiz isso antes de criar a classe Curve2D e fiquei com preguiça de 
## mudar o código para usar ela.
class Questao3(Scene):

    def construct(self):
        # Two elipses, parametrized in [-2pi, 2pi]
        a = 2
        b = 1
        init_t = -2 * PI

        # Time parameter
        time = ValueTracker(init_t)
        time.add_updater(lambda mobject, dt: mobject.increment_value(dt))

        # Symbolic time parameter
        t = symbols('t')

        speed = .5

        ## Anti clockwise parametrization

        name_ac = MathTex(r"\alpha(t) = (a \cos(t), b \sin(t))").set_x(-3).set_y(3)

        # Symbolic coordinates
        x_sym_ac = a*cos(t*speed) - 3
        y_sym_ac = b*sin(t*speed)

        # Anti clockwise curve
        x_ac = lambdify(t, x_sym_ac, "numpy")
        y_ac = lambdify(t, y_sym_ac, "numpy")
        x_prime_ac = lambdify(t, diff(x_sym_ac, t), "numpy")
        y_prime_ac = lambdify(t, diff(y_sym_ac, t), "numpy")


        def parametrization_ac(t):
            return np.array([x_ac(t), y_ac(t), 0])

        # Starting location
        start_ac = parametrization_ac(init_t)

        def velocity_ac(t):
            return np.array([x_prime_ac(t), y_prime_ac(t), 0])

        def curvature_ac(t):
            return a*b / (a**2 * np.sin(t)**2 + b**2 * np.cos(t)**2)**(3/2)


        dot_ac = Dot(point=start_ac, color=RED_C)
        dot_ac.add_updater(lambda x: x.move_to(parametrization_ac(time.get_value())))

        vec_ac = Vector(velocity_ac(0), color=RED_C,
                        max_tip_length_to_length_ratio=100)
        vec_ac.add_updater(lambda x: x.put_start_and_end_on(dot_ac.get_center(),
                                                            velocity_ac(time.get_value()) + dot_ac.get_center()))

        path_ac = TracedPath(dot_ac.get_center, stroke_width=3, stroke_color=RED_C)

        # Curvature section

        on_screen_curv_ac = Variable(a/b**2, MathTex("k_{\\alpha}(t)"),
                                     num_decimal_places=2).scale(.7)
        def curv_updater_ac(m):
            m.move_to(parametrization_ac(time.get_value()) + velocity_ac(time.get_value())*1.5)
            m.tracker.set_value(curvature_ac(time.get_value()))

        on_screen_curv_ac.add_updater(curv_updater_ac)

        ## Clockwise parametrization

        name_c = MathTex(r"\beta(t) = \alpha(-t)").set_x(3).set_y(3)

        # Symbolic coordinates
        x_sym_c = a*cos(-t*speed) + 3
        y_sym_c = b*sin(-t*speed)

        # Clockwise curve
        x_c = lambdify(t, x_sym_c, "numpy")
        y_c = lambdify(t, y_sym_c, "numpy")
        x_prime_c = lambdify(t, diff(x_sym_c, t), "numpy")
        y_prime_c = lambdify(t, diff(y_sym_c, t), "numpy")


        def parametrization_c(t):
            return np.array([x_c(t), y_c(t), 0])

        # Starting location
        start_c = parametrization_c(init_t)

        def velocity_c(t):
            return np.array([x_prime_c(t), y_prime_c(t), 0])

        def curvature_c(t):
            return -curvature_ac(-t)

        dot_c = Dot(point=start_c, color=RED_C)
        dot_c.add_updater(lambda x: x.move_to(parametrization_c(time.get_value())))

        vec_c = Vector(velocity_c(0), color=RED_C,
                     max_tip_length_to_length_ratio=100)
        vec_c.add_updater(lambda x: x.put_start_and_end_on(dot_c.get_center(),
                                                         velocity_c(time.get_value()) + dot_c.get_center()))

        path_c = TracedPath(dot_c.get_center, stroke_width=3, stroke_color=RED_C)

        # Curvature section

        on_screen_curv_c = Variable(a/b**2, MathTex("k_{\\beta}(t)"),
                                     num_decimal_places=2).scale(.7)
        def curv_updater_c(m):
            m.move_to(parametrization_c(time.get_value()) + velocity_c(time.get_value())*1.5)
            m.tracker.set_value(curvature_c(time.get_value()))

        on_screen_curv_c.add_updater(curv_updater_c)

        self.add(time,
                 name_ac, dot_ac, path_ac, vec_ac, on_screen_curv_ac,
                 name_c, dot_c, path_c, vec_c, on_screen_curv_c)
        self.wait(4 * PI)


class Questao4(Scene):

    def construct(self):

        speed = .2

        # Regular curve
        init_t = -4

        t = Symbol("t")
        x = t
        y = cosh(t)
        curve = Curve2D("alpha", t, [x, y], init_t, RED_C, speed=.5, scale=1,
                        offset=[-3, -4])

        # Unit speed
        x_us = asinh(t)
        y_us = cosh(asinh(t))
        curve_us = Curve2D("beta", t, [x_us, y_us], init_t, BLUE_C, speed=.5,
                           scale=1, offset=[3, -4])

        title = Text("Catenária").set_y(-3).scale(.5)

        self.add(title,
                 curve.time, curve.dot_manim, curve.path_manim,
                 curve.velocity_manim, curve.acceleration_manim,
                 curve.velocity_name_manim, curve.acceleration_name_manim,
                 #
                 curve_us.time, curve_us.dot_manim, curve_us.path_manim,
                 curve_us.velocity_manim, curve_us.acceleration_manim,
                 curve_us.velocity_name_manim, curve_us.acceleration_name_manim)

        self.wait(8)


class Questao5(Scene):

    def construct(self):

        init_t = -6
        t = Symbol("t")
        x = t
        y = sin(t)**3

        circle = manim.Circle(radius=1, color=BLUE_C)

        curve = Curve2D("alpha", t, [x, y], init_t, RED_C, speed=1, scale=1,
                        offset=[0, 0])

        # Tangent Indicatrix
        vel = Vector(curve.velocity(init_t), color=BLUE_C)
        def vel_updater(m):
            m.put_start_and_end_on([0, 0, 0], curve.velocity(curve.time.get_value()))
        vel.add_updater(vel_updater)

        def tangent_indicatrix(t):
            v = curve.velocity(t)
            return v / np.linalg.norm(v)

        tan_ind = Dot(tangent_indicatrix(init_t))
        tan_ind.add_updater(lambda x:
                            x.move_to(tangent_indicatrix(curve.time.get_value())))

        line = manim.Line([1, 10, 0], [1, -10, 0])

        title = Text(r"Tangente Indicatrix").scale(.5).set_y(3.2)

        self.add(NumberPlane(), vel, circle, tan_ind, line,
                 curve.time, curve.dot_manim, curve.path_manim,
                 curve.velocity_manim, curve.velocity_name_manim, title)

        self.wait(12)


class Questao10(Scene):


    def construct(self):

        def func(t):
            return np.array([np.sin(2 * t), np.sin(3 * t), 0])

        curve = ParametricFunction(func, t_range=np.array([0, TAU]), fill_opacity=0).set_color(RED).scale(2)
        self.wait()
        self.play(Create(curve))
        self.play(curve.animate.rotate(PI/4))
        self.play(curve.animate.shift(2 * LEFT))
        self.play(curve.animate.shift(2 * RIGHT))
        self.play(curve.animate.rotate(PI/4))
        self.play(curve.animate.shift(2 * RIGHT + 2 * UP))
        self.play(curve.animate.shift(2 * LEFT + 2 * DOWN))
        self.play(curve.animate.rotate(PI/2))
        self.wait()

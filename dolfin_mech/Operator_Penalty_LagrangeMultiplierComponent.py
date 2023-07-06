#coding=utf8

################################################################################
###                                                                          ###
### Created by Martin Genet, 2018-2023                                       ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
<<<<<<< HEAD:dolfin_mech/Operator_TensorSymmetry.py
###                                                                          ###
### And Mahdi Manoochehrtayebi, 2021-2023                                    ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
=======
>>>>>>> f86a0353825e9b67fbb7f784495eb1ae7db1ae7b:dolfin_mech/Operator_Penalty_LagrangeMultiplierComponent.py
################################################################################

import dolfin

import dolfin_mech as dmech
from .Operator import Operator

################################################################################

class LagrangeMultiplierComponentPenaltyOperator(Operator):

    def __init__(self,
            lambda_bar,
            lambda_bar_test,
            i, j,
            measure,
            pen_val=None, pen_ini=None, pen_fin=None):

        self.measure = measure

        self.tv_pen = dmech.TimeVaryingConstant(
            val=pen_val, val_ini=pen_ini, val_fin=pen_fin)
        pen = self.tv_pen.val

        self.res_form = pen * lambda_bar[i,j] * lambda_bar_test[i,j] * self.measure



    def set_value_at_t_step(self,
            t_step):

        self.tv_pen.set_value_at_t_step(t_step)

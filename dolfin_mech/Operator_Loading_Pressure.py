#coding=utf8

################################################################################
###                                                                          ###
### Created by Martin Genet, 2018-2022                                       ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import dolfin

import dolfin_mech as dmech
from .Operator import Operator

################################################################################

class PressureLoadingOperator(Operator):

    def __init__(self,
            U_test,
            kinematics,
            N,
            measure,
            P=None,
            P_ini=None,
            P_fin=None):

        self.measure = measure

        self.tv_P = dmech.TimeVaryingConstant(
            val=P, val_ini=P_ini, val_fin=P_fin)
        P = self.tv_P.val

        T = dolfin.dot(-P * N, dolfin.inv(kinematics.F))
        self.res_form = - dolfin.inner(T, U_test) * kinematics.J * self.measure



    def set_value_at_t_step(self,
            t_step):

        self.tv_P.set_value_at_t_step(t_step)

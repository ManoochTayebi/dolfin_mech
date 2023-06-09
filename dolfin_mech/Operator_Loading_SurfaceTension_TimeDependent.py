#coding=utf8

################################################################################
###                                                                          ###
### Created by Martin Genet, 2018-2022                                       ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import dolfin
import math

import dolfin_mech as dmech
from .Operator import Operator

################################################################################

class SurfaceTensionLoadingOperatorTimeDependent(Operator):

    def __init__(self,
            C,
            k1,
            k2,
            gamma_0,
            Gamma_hat,
            m1,
            Gamma_history,
            qois,
            U,
            U_test,
            kinematics,
            N,
            measure,
            dS,
            gamma_val=None, gamma_ini=None, gamma_fin=None):

        self.measure = measure
        self.N = N
        self.dS = dS
        self.Gamma_history = Gamma_history
        self.C = C
        self.k1 = k1
        self.k2 = k2
        self.gamma_0 = gamma_0
        self.Gamma_hat = Gamma_hat
        self.m1 = m1

        self.tv_gamma = dmech.TimeVaryingConstant(
            val=gamma_val, val_ini=gamma_ini, val_fin=gamma_fin)
        gamma = self.tv_gamma.val

        self.tv_dt = dmech.TimeVaryingConstant(0.)
        dt = self.tv_dt.val

        self.Gamma_dt = dmech.TimeVaryingConstant(0.)
        Gamma = self.Gamma_dt.val

        self.gamma_t_val = dmech.TimeVaryingConstant(0.)
        gamma_t = self.gamma_t_val.val


        self.S_t = dmech.TimeVaryingConstant(0.)
        S = self.S_t.val


        # # self.tv_gamma.surface_change_rate(kinematics, dt)
        # # print("gamma =" +str(gamma))
        # FmTN = dolfin.dot(dolfin.inv(kinematics.F).T, N)
        # T = dolfin.sqrt(dolfin.inner(FmTN, FmTN))
        # # Pi = gamma * T * kinematics.J * self.measure
        # Pi = gamma_t * T * kinematics.J * self.measure
        # self.res_form = dolfin.derivative(Pi, U, U_test)

        dim = U.ufl_shape[0]
        I = dolfin.Identity(dim)
        Pi = gamma_t * (1 + dolfin.inner(
            kinematics.E,
            I - dolfin.outer(N,N))) * self.measure
        self.res_form = dolfin.derivative(Pi, U, U_test) # MG20211220: Is that correct?!

        self.kinematics=kinematics

    # def surface_change_rate(self,
    #         dt):

    #     self.tv_gamma.surface_change_rate(self.kinematics, dt)

    def set_value_at_t_step(self,
            t_step):

        self.tv_gamma.set_value_at_t_step(t_step)
        print("t_step =" +str(t_step))
        # print("value at t_step = " +str(self.tv_gamma.set_value_at_t_step(t_step)))

    def returne_surface_rate(self):
        self.tv_gamma.surface_change_rate()


    def set_dt(self,
            t_step,
            dt):
        S_old = dolfin.assemble(dolfin.sqrt(dolfin.dot(dolfin.inv(self.kinematics.F_old.T), self.N)[0]**2+ dolfin.dot(dolfin.inv(self.kinematics.F_old.T), self.N)[1]**2)*self.kinematics.J_old * self.dS(0))
        S = dolfin.assemble(dolfin.sqrt(dolfin.dot(dolfin.inv(self.kinematics.F.T), self.N)[0]**2+ dolfin.dot(dolfin.inv(self.kinematics.F.T), self.N)[1]**2)*self.kinematics.J * self.dS(0))
        dS_dt = (S - S_old)/dt
        print("S_old: " +str(S_old))
        print("S: " +str(S))
        print("dS_dt: "+str(dS_dt))
        print("dt: "+str(dt))
        self.tv_dt.set_value(dt)
        self.Gamma_dt.set_value(dS_dt)
        self.S_t.set_value(S)
        

        t0 = t_step - dt
        a = self.k1*self.C + self.k2/S + dS_dt/S
        b = self.k1 * self.C * self.Gamma_hat
        # A = (self.Gamma_history[-1] - b/a) * math.e**(-a*t0)
        A = (self.Gamma_history[-1] - b/a)

        Gamma = b/a + A * math.e**(-a*dt)
        self.Gamma_history.append(Gamma)

        gamma_new = self.gamma_0 - self.m1 * Gamma/self.Gamma_hat
        print("gamma: " +str(gamma_new))
        self.gamma_t_val.set_value(gamma_new)

################################################################################

class SurfaceTension0LoadingOperator(Operator):

    def __init__(self,
            u,
            u_test,
            kinematics,
            N,
            measure,
            gamma_val=None, gamma_ini=None, gamma_fin=None):

        self.measure = measure

        self.tv_gamma = dmech.TimeVaryingConstant(
            val=gamma_val, val_ini=gamma_ini, val_fin=gamma_fin)
        gamma = self.tv_gamma.val

        dim = u.ufl_shape[0]
        I = dolfin.Identity(dim)
        Pi = gamma * (1 + dolfin.inner(
            kinematics.E,
            I - dolfin.outer(N,N))) * self.measure
        self.res_form = dolfin.derivative(Pi, u, u_test) # MG20211220: Is that correct?!



    def set_value_at_t_step(self,
            t_step):

        self.tv_gamma.set_value_at_t_step(t_step)

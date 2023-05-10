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



class SurfacePressureGradientLoadingOperator(Operator):

    def __init__(self,
            X,
            x0,
            x0_test,
            kinematics,
            U,
            U_test,
            Phis0,
            V0, 
            N,
            lbda,
            lbda_test,
            mu,
            mu_test,
            p,
            p_test,
            gamma,
            gamma_test,
            dS,
            dV,
            rho_solid,
            proj_op,
            # Xm,
            P0_val=None, P0_ini=None, P0_fin=None,
            f_val=None, f_ini=None, f_fin=None):


        self.measure = dV
        self.dS = dS

        self.tv_F = dmech.TimeVaryingConstant(
            val=f_val, val_ini=f_ini, val_fin=f_fin)
        F = self.tv_F.val

        self.tv_P0 = dmech.TimeVaryingConstant(
            val=P0_val, val_ini=P0_ini, val_fin=P0_fin)
        P0 = self.tv_P0.val

        self.tv_fy = dmech.TimeVaryingConstant(
            val=None, val_ini=0., val_fin=f_fin[1])
        fy = self.tv_fy.val


        nf = dolfin.dot(N, dolfin.inv(kinematics.F))
        nf_norm = (1/2 * dolfin.inner(nf,nf))**(1/2)
        n = nf/nf_norm

        self.V0 = dolfin.Constant(V0)
        rho0 = dolfin.assemble(Phis0 * self.measure)

        x = X + U
        x_tilde = x-x0

        

        Pconst = P0 - rho_solid * fy * ( x[1]- x0[1])

        # grads_p = dolfin.dot(proj_op, dolfin.dot(dolfin.grad(p-P0), proj_op))
        # grads_p_test = dolfin.dot(proj_op, dolfin.dot(dolfin.grad(p_test), proj_op))
        # grads_p = dolfin.dot(proj_op, dolfin.dot(dolfin.inv(kinematics.F).T * dolfin.grad(p-P0), proj_op))
        # grads_p_test = dolfin.dot(proj_op, dolfin.dot(dolfin.inv(kinematics.F).T * dolfin.grad(p_test), proj_op))
        grads_p = dolfin.dot(dolfin.grad(p-Pconst), dolfin.inv(kinematics.F)) - n*(dolfin.dot(n,dolfin.dot(dolfin.grad(p-Pconst), dolfin.inv(kinematics.F))))
        grads_p_test = dolfin.dot(dolfin.grad(p_test), dolfin.inv(kinematics.F)) - n*(dolfin.dot(n,dolfin.dot(dolfin.grad(p_test), dolfin.inv(kinematics.F))))


        self.res_form = dolfin.Constant(1e-8)*p*p_test * kinematics.J * self.measure
        self.res_form -= dolfin.inner(rho_solid*Phis0*F, U_test) * self.measure
        self.res_form -=  dolfin.inner(-p * n, U_test) * nf_norm * kinematics.J * dS
        self.res_form += dolfin.inner(rho_solid*Phis0*F, lbda_test) * self.measure
        self.res_form += dolfin.inner(-p * n, lbda_test) * nf_norm * kinematics.J * dS
        self.res_form += -dolfin.dot(lbda, n) * p_test * nf_norm * kinematics.J * self.dS
        self.res_form += - dolfin.dot(mu, dolfin.cross(x_tilde, n)) *  p_test * nf_norm * kinematics.J * self.dS
        self.res_form += gamma  *  p_test * nf_norm * kinematics.J * self.dS
        self.res_form +=  dolfin.inner(grads_p, grads_p_test) * nf_norm * kinematics.J * self.dS
        self.res_form += dolfin.inner(dolfin.cross(x_tilde, -p * n), mu_test) * nf_norm * kinematics.J * dS
        self.res_form += (p-Pconst)*gamma_test * nf_norm * kinematics.J * dS
        self.res_form -= dolfin.inner((Phis0 * x / dolfin.Constant(rho0) - x0/self.V0), x0_test) * self.measure

         

    def set_value_at_t_step(self,
            t_step):       
        self.tv_F.set_value_at_t_step(t_step)
        self.tv_P0.set_value_at_t_step(t_step)
        self.tv_fy.set_value_at_t_step(t_step)


################################################################################

class SurfacePressureGradient0LoadingOperator(Operator):

    def __init__(self,
            x,
            x0,
            n,
            u_test,
            lbda,
            lbda_test,
            mu,
            mu_test,
            p,
            p_test,
            gamma,
            gamma_test,
            dS,
            dV,
            rho_solid,
            phis,
            proj_op,
            # Xm,
            P0_val=None, P0_ini=None, P0_fin=None,
            f_val=None, f_ini=None, f_fin=None):


        self.dS = dS
        self.measure = dV

        self.tv_f = dmech.TimeVaryingConstant(
            val=f_val, val_ini=f_ini, val_fin=f_fin)
        f = self.tv_f.val

        self.tv_P0 = dmech.TimeVaryingConstant(
            val=P0_val, val_ini=P0_ini, val_fin=P0_fin)
        P0 = self.tv_P0.val

        self.tv_fy = dmech.TimeVaryingConstant(
            val=None, val_ini=0, val_fin=f_fin[1])
        fy = self.tv_fy.val

        x_tilde = x-dolfin.Constant(x0)


        Pconst = P0 - rho_solid*fy * ( x[1]- dolfin.Constant(x0[1]))
        

        # p = P0 + dolfin.dot(lbda, n) + dolfin.dot(mu, dolfin.cross(x_tilde, n))

        # self.res_form = - dolfin.inner(rho_solid*phis*f, u_test) * self.measure
        # self.res_form -= dolfin.inner(-P * n, u_test) * self.dS
        # self.res_form += dolfin.inner(dolfin.cross(x_tilde, -P * n), mu_test) * self.dS
        # self.res_form += dolfin.inner(rho_solid*phis*f, lbda_test) * self.measure
        # self.res_form += dolfin.inner(-P * n, lbda_test) * self.dS



        # grads_p = dolfin.dot(proj_op, dolfin.dot(dolfin.grad(p-P0), proj_op))
        # grads_p_test = dolfin.dot(proj_op, dolfin.dot(dolfin.grad(p_test), proj_op))
        grads_p = dolfin.grad(p-Pconst) - n*(dolfin.dot(n,dolfin.grad(p-Pconst)))
        grads_p_test = dolfin.grad(p_test) - n*(dolfin.dot(n,dolfin.grad(p_test)))


        self.res_form = dolfin.Constant(1e-8)*p*p_test * self.measure
        self.res_form -= dolfin.inner(rho_solid*phis*f, u_test) * self.measure
        self.res_form -=  dolfin.inner(-p * n, u_test) * dS
        self.res_form += dolfin.inner(rho_solid*phis*f, lbda_test) * self.measure
        self.res_form += dolfin.inner(-p * n, lbda_test) * dS
        self.res_form += -dolfin.dot(lbda, n) *  p_test * self.dS
        self.res_form += - dolfin.dot(mu, dolfin.cross(x_tilde, n)) *  p_test * self.dS
        self.res_form += gamma  *  p_test * self.dS
        self.res_form +=  dolfin.inner(grads_p, grads_p_test) * self.dS
        self.res_form += dolfin.inner(dolfin.cross(x_tilde, -p * n), mu_test) * dS
        self.res_form += (p-Pconst)*gamma_test * dS



    def set_value_at_t_step(self,
            t_step):
        self.tv_f.set_value_at_t_step(t_step)
        self.tv_P0.set_value_at_t_step(t_step)
        self.tv_fy.set_value_at_t_step(t_step)
        

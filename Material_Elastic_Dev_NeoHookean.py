#coding=utf8

################################################################################
###                                                                          ###
### Created by Martin Genet, 2018-2020                                       ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

# from builtins import *

import dolfin

import dolfin_cm as dcm
from .Material_Elastic_Dev import DevElasticMaterial

################################################################################

class NeoHookeanDevElasticMaterial(DevElasticMaterial):



    def __init__(self,
            parameters):

        if ("C1" in parameters) and ("C2" in parameters):
            self.C1 = dolfin.Constant(parameters["C1"])
        elif ("mu" in parameters):
            self.mu = dolfin.Constant(parameters["mu"])
            self.C1 = self.mu/2
        elif ("E" in parameters) and ("nu" in parameters):
            self.E  = dolfin.Constant(parameters["E"])
            self.nu = dolfin.Constant(parameters["nu"])
            self.mu = self.E/2/(1+self.nu)
            self.C1 = self.mu/2



    def get_free_energy(self,
            U=None,
            C=None):

        if (C is None):
            dim = U.ufl_shape[0]
            I = dolfin.Identity(dim)
            F = I + dolfin.grad(U)
            C = F.T * F
        else:
            assert (C.ufl_shape[0] == C.ufl_shape[1])
            dim = C.ufl_shape[0]
            I = dolfin.Identity(dim)

        JF    = dolfin.sqrt(dolfin.det(C))
        IC    = dolfin.tr(C)
        C_inv = dolfin.inv(C)

        if   (dim == 2):
            Psi   =   self.C1 * (IC - 2 - 2*dolfin.ln(JF)) #MG20200206: plane strain
            Sigma = 2*self.C1 * (I - C_inv)                #MG20200206: plane strain
        elif (dim == 3):
            Psi   =   self.C1 * (IC - 3 - 2*dolfin.ln(JF))
            Sigma = 2*self.C1 * (I - C_inv)

        return Psi, Sigma

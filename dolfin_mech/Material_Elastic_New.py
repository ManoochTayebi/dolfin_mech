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
from .Material_Elastic import ElasticMaterial

################################################################################

class NewElasticMaterial(ElasticMaterial):



    def __init__(self,
            kinematics,
            parameters,
            decoup=False):

        self.kinematics = kinematics

        self.C1 = dolfin.Constant(parameters["C1"])
        self.C2 = dolfin.Constant(parameters["C2"])

        if   (self.kinematics.dim == 2):
                self.Psi   = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)) + 300*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                self.Sigma = self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv*(1 + self.kinematics.IC)) * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)) + 2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # Mahdi

                self.Psi_old   = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old - 3))) + 300*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                self.Sigma_old = self.C1 * self.kinematics.J_old**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv_old*(1 + self.kinematics.IC_old)) * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC_old - 3))) + 2*300*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old # Mahdi


                # assert (0), "ToDo. Aborting."
        elif (self.kinematics.dim == 3):
            self.Psi = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.IC - 3)) + 300*self.C1 *(self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
            self.Sigma = self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3 *self.kinematics.IC*self.kinematics.C_inv) * dolfin.exp(self.C2*(self.C2(self.kinematics.IC - 3))) + 2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv

            self.Psi_old = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.IC_old - 3)) + 300*self.C1 *(self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
            self.Sigma_old = self.C1 * self.kinematics.J_old**(-2/3) *(self.kinematics.I - 1/3 *self.kinematics.IC_old*self.kinematics.C_inv_old) * dolfin.exp(self.C2*(self.C2(self.kinematics.IC_old - 3))) + 2*300*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old




        self.P = self.kinematics.F * self.Sigma
        self.P_old = self.kinematics.F_old * self.Sigma_old # Mahdi 

        self.sigma = self.P * self.kinematics.F.T / self.kinematics.J
        self.sigma_old = self.P_old * self.kinematics.F_old.T / self.kinematics.J_old # Mahdi
        
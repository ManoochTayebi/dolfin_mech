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

class CiarletGeymonatElasticMaterial(ElasticMaterial):



    def __init__(self,
            kinematics,
            parameters,
            decoup=False):

        self.kinematics = kinematics

        self.lmbda = self.get_lambda_from_parameters(parameters)
        self.mu = self.get_mu_from_parameters(parameters)

        self.K = (3*self.lmbda + 2*self.mu)/3

        # if ("beta" in parameters) and ("delta" in parameters):
        #     self.beta = dolfin.Constant(parameters["beta"])
        #     self.delta = dolfin.Constant(parameters["delta"])
        # else:
        #     assert (0),\
        #         "No parameter found: \"+str(parameters)+\". Must provide beta & delta. Aborting."

        # self.Psi = (self.beta) * (dolfin.exp(self.delta*(self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))) - 1)

        # self.Sigma = (self.beta) * dolfin.exp(self.delta*(self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))) * (2*self.delta) * (self.kinematics.J**2 - 1) * self.kinematics.C_inv


        # self.Sigma_old = (self.beta) * dolfin.exp(self.delta*(self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J))) * (2*self.delta) * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old

        

        # # self.P = dolfin.diff(self.Psi, self.kinematics.F) # MG20220426: Cannot do that for micromechanics problems
        # self.P = self.kinematics.F * self.Sigma
        # self.sigma = self.P * self.kinematics.F.T / self.kinematics.J

        # self.P_old = self.kinematics.F_old * self.Sigma_old
        # self.sigma_old = self.P_old * self.kinematics.F_old.T / self.kinematics.J_old
        # self.Sigma_zz = 1

        if (decoup):
            # self.Psi = (self.lmbda/4) * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J)) # MG20180516: In 2d, plane strain

            # self.Sigma = (self.lmbda/2) * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # MG20200206: Cannot differentiate Psi wrt to C because J is not defined as a function of C
            # self.Sigma_old = (self.lmbda/2) * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old # Mahdi


            self.Psi = (self.lmbda/4 + self.mu/6) * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J)) # MG20180516: In 2d, plane strain

            self.Sigma = (self.lmbda/2 + self.mu/3) * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # MG20200206: Cannot differentiate Psi wrt to C because J is not defined as a function of C
            self.Sigma_old = (self.lmbda/2 + self.mu/3) * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old # Mahdi

            self.Sigma_zz = (self.lmbda/2 + self.mu/3) * (self.kinematics.J**2 - 1)

        else:
            self.Psi = (self.lmbda/4) * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J)) # MG20180516: In 2d, plane strain
  
            self.Sigma = (self.lmbda/2) * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # MG20200206: Cannot differentiate Psi wrt to C because J is not defined as a function of C
            self.Sigma_old = (self.lmbda/2) * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old # Mahdi

            self.Sigma_zz = (self.lmbda/2) * (self.kinematics.J**2 - 1)


        # self.P = dolfin.diff(self.Psi, self.kinematics.F) # MG20220426: Cannot do that for micromechanics problems
        # self.P = (self.lmbda/2) * (self.kinematics.J**2 - 1) * self.kinematics.F_inv.T
        self.P = self.kinematics.F * self.Sigma
        self.P_old = self.kinematics.F_old * self.Sigma_old # Mahdi

        self.sigma = self.P * self.kinematics.F.T / self.kinematics.J
        self.sigma_old = self.P_old * self.kinematics.F_old.T / self.kinematics.J_old # Mahdi
        

    # def get_free_energy(self,
    #         U=None,
    #         C=None):

    #     C  = self.get_C_from_U_or_C(U, C)
    #     JF = dolfin.sqrt(dolfin.det(C)) # MG20200207: Watch out! This is well defined for inverted elements!

    #     Psi   = (self.lmbda/4) * (JF**2 - 1 - 2*dolfin.ln(JF)) # MG20180516: in 2d, plane strain
    #     Sigma = 2*dolfin.diff(Psi, C)

    #     # C_inv = dolfin.inv(C)
    #     # Sigma = (self.lmbda/2) * (JF**2 - 1) * C_inv

    #     return Psi, Sigma

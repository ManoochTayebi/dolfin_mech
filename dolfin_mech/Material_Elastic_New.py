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
        self.C3 = dolfin.Constant(parameters["C3"])


        if   (self.kinematics.dim == 2):
                # self.Psi   = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3))\
                #            + 300*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                
                # self.Sigma = self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv*(1 + self.kinematics.IC))* dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3))\
                #            + 2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # Mahdi



                # self.Psi_old   = self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3))\
                #            + 300*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                
                # self.Sigma_old = self.C1 * self.kinematics.J_old**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv_old*(1 + self.kinematics.IC_old)) * dolfin.exp(self.C2*(self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3))\
                #            + 2*300*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old # Mahdi



                # self.Psi   = self.C1/self.C2/4 * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**2)\
                #            +  self.C1/10 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            + 300*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                
                # self.Sigma = self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv*(1 + self.kinematics.IC)) *(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)* dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**2)\
                #            +  2 * self.C1/10 * (self.kinematics.I - self.kinematics.C_inv)\
                #            + 2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # Mahdi



                # self.Psi_old   = self.C1/self.C2/4 * dolfin.exp(self.C2*(self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3)**2)\
                #            +  self.C1/10 * (self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))\
                #            + 300*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                
                # self.Sigma_old= self.C1 * self.kinematics.J_old**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv_old*(1 + self.kinematics.IC_old)) *(self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3)* dolfin.exp(self.C2*(self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3)**2)\
                #            +  2 * self.C1/10 * (self.kinematics.I - self.kinematics.C_inv_old)\
                #            + 2*300*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old # Mahdi

                # self.Psi   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                #            +  self.C1/10 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv) * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                #            +  2 * self.C1/10 * (self.kinematics.I - self.kinematics.C_inv)\
                #            +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv

                # self.Psi_old   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**self.C3) \
                #                +  self.C1/10 * (self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))\
                #                +  100*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                # self.Sigma_old =  self.C1 * (self.kinematics.I - self.kinematics.C_inv_old) * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**self.C3 ) \
                #                +  2 * self.C1/10 * (self.kinematics.I - self.kinematics.C_inv_old)\
                #                +  2*100*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old



                # self.Psi   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*((self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**1.2)) +  300*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv*(1 + self.kinematics.IC)) * dolfin.exp(self.C2*((self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**1.2)) +  2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv

                # self.Psi_old   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*((self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3)**1.2)) +  300*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                # self.Sigma_old =  self.C1 * self.kinematics.J_old**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv_old*(1 + self.kinematics.IC_old)) * dolfin.exp(self.C2*((self.kinematics.J_old**(-2/3) * (1 + self.kinematics.IC_old) - 3)**1.2)) +  2*300*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old




                # self.Psi   =  self.C1/self.C2/4 * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                #            +  self.C1/1 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv)* (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J)) * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                #            +  2 * self.C1/1 * (self.kinematics.I - self.kinematics.C_inv)\
                #            +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv 

                # self.Psi_old   =  self.C1/self.C2/4 * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**self.C3) \
                #                +  self.C1/1 * (self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))\
                #                +  100*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                # self.Sigma_old =  self.C1 * (self.kinematics.I - self.kinematics.C_inv_old)*self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old)) * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**self.C3) \
                #                +  2 * self.C1/1 * (self.kinematics.I - self.kinematics.C_inv_old)\
                #                +  2*100*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old


                ############ Power 2

                # self.Psi   =  self.C1/self.C2/4 * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**2) \
                #            +  self.C1 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv) * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J)) * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**2) \
                #            +  2 * self.C1 * (self.kinematics.I - self.kinematics.C_inv)\
                #            +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv 
                

                # self.Psi   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                #            +  self.C1/1 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv)* (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**(self.C3 - 1) * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                #            +  2 * self.C1/1 * (self.kinematics.I - self.kinematics.C_inv)\
                #            +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv 

                # self.Psi_old   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**self.C3) \
                #                +  self.C1/1 * (self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))\
                #                +  100*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                # self.Sigma_old =  self.C1 * (self.kinematics.I - self.kinematics.C_inv_old)*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**(self.C3 - 1) * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**self.C3) \
                #                +  2 * self.C1/1 * (self.kinematics.I - self.kinematics.C_inv_old)\
                #                +  2*100*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old


                ############ Power 1.5

                # self.Psi   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*dolfin.sqrt(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**2) \
                #            +  self.C1/10 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv) * dolfin.exp(self.C2*dolfin.sqrt(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**2) \
                #            +  2 * self.C1/10 * (self.kinematics.I - self.kinematics.C_inv)\
                #            +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv 

                # self.Psi_old   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*dolfin.sqrt(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**2) \
                #                +  self.C1/10 * (self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))\
                #                +  100*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                # self.Sigma_old =  self.C1 * (self.kinematics.I - self.kinematics.C_inv_old) * dolfin.exp(self.C2*dolfin.sqrt(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))**2) \
                #                +  2 * self.C1/10 * (self.kinematics.I - self.kinematics.C_inv_old)\
                #                +  2*100*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old
                

                ############ Power 1


                # self.Psi   =  self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))) \
                #            +  self.C1 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                # self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv) * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))) \
                #            +  2 * self.C1 * (self.kinematics.I - self.kinematics.C_inv)\
                #            +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv 

                # self.Psi_old   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))) \
                #                +  self.C1/1 * (self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))\
                #                +  100*self.C1 * (self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
                # self.Sigma_old =  self.C1 * (self.kinematics.I - self.kinematics.C_inv_old) * dolfin.exp(self.C2*(self.kinematics.IC_old - 2 - 2*dolfin.ln(self.kinematics.J_old))) \
                #                +  2 * self.C1/1 * (self.kinematics.I - self.kinematics.C_inv_old)\
                #                +  2*100*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old


                ################# IC_bar


                # self.Psi   = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3))\
                #            + 300*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                
                # self.Sigma = self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv*(1 + self.kinematics.IC))* dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3))\
                #            + 2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # Mahdi
                
                # self.Psi   = self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**self.C3)\
                #            +  self.C1/1 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                #            + 100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                
                # self.Sigma = self.C1 * (self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**(self.C3 - 1) * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3*self.kinematics.C_inv*(1 + self.kinematics.IC)) * dolfin.exp(self.C2*(self.kinematics.J**(-2/3) * (1 + self.kinematics.IC) - 3)**self.C3)\
                #            +  2 * self.C1/1 * (self.kinematics.I - self.kinematics.C_inv)\
                #            + 2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv # Mahdi
                
                
                ################### Power C3

                self.Psi   =  self.C1/self.C2/self.C3/2 * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                           +  8*self.C1 * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))\
                           +  100*self.C1 * (self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
                self.Sigma =  self.C1 * (self.kinematics.I - self.kinematics.C_inv) * (self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**(self.C3 - 1) * dolfin.exp(self.C2*(self.kinematics.IC - 2 - 2*dolfin.ln(self.kinematics.J))**self.C3) \
                           +  2 * 8*self.C1 * (self.kinematics.I - self.kinematics.C_inv)\
                           +  2*100*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv 

                # assert (0), "ToDo. Aborting."
        elif (self.kinematics.dim == 3):
            self.Psi = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.IC - 3)) + 300*self.C1 *(self.kinematics.J**2 - 1 - 2*dolfin.ln(self.kinematics.J))
            self.Sigma = self.C1 * self.kinematics.J**(-2/3) *(self.kinematics.I - 1/3 *self.kinematics.IC*self.kinematics.C_inv) * dolfin.exp(self.C2*(self.C2(self.kinematics.IC - 3))) + 2*300*self.C1 * (self.kinematics.J**2 - 1) * self.kinematics.C_inv

            self.Psi_old = self.C1/self.C2/2 * dolfin.exp(self.C2*(self.kinematics.IC_old - 3)) + 300*self.C1 *(self.kinematics.J_old**2 - 1 - 2*dolfin.ln(self.kinematics.J_old))
            self.Sigma_old = self.C1 * self.kinematics.J_old**(-2/3) *(self.kinematics.I - 1/3 *self.kinematics.IC_old*self.kinematics.C_inv_old) * dolfin.exp(self.C2*(self.C2(self.kinematics.IC_old - 3))) + 2*300*self.C1 * (self.kinematics.J_old**2 - 1) * self.kinematics.C_inv_old




        self.P = self.kinematics.F * self.Sigma
        # self.P_old = self.kinematics.F_old * self.Sigma_old # Mahdi 

        self.sigma = self.P * self.kinematics.F.T / self.kinematics.J
        # self.sigma_old = self.P_old * self.kinematics.F_old.T / self.kinematics.J_old # Mahdi
        
import numpy as np
import json

class BasicInfo():
    def __init__(self, Nx, Ny, dx, dy, Nt, dt):
        self.Nx = Nx
        self.Ny = Ny
        self.dx = dx
        self.dy = dy
        self.Nt = Nt
        self.dt = dt
        self.t = np.linspace(0, (Nt-1)*dt, Nt)

class Source(BasicInfo):

    def __init__(self, Nx, Ny, dx, dy, Nt, dt):
        self.source_num = 0
        self.source_position = []
        self.source_fn = []
        BasicInfo.__init__(self, Nx, Ny, dx, dy, Nt, dt)

    def SetSource(self, source_num, source_position, fre, center):
        self.source_num = source_num
        self.source_position = np.array(source_position)
        self.source_fn = np.zeros([self.source_num, self.Nt])
        for i in range(self.source_num):
            self.source_fn[i,:] = self.RickerFn(fre, center)

    def RickerFn(self, fre, center):
        return (1 - 2 * np.pi**2 * fre**2 * (self.t-center)**2) * np.exp(-1 * np.pi**2 * fre**2 * (self.t-center)**2)
    
class Receiver(BasicInfo):
    
    def __init__(self, Nx, Ny, dx, dy, Nt, dt):
        self.receiver_num = 0
        self.receiver_position = []
        BasicInfo.__init__(self, Nx, Ny, dx, dy, Nt, dt)

    def SetReceiver(self, receiver_num, receiver_position):
        self.receiver_num = receiver_num
        self.receiver_position = np.array(receiver_position)

class Model(Source, Receiver, BasicInfo):

    def __init__(self, Nx, Ny, dx, dy, Nt, dt):
        Source.__init__(self, Nx, Ny, dx, dy, Nt, dt)
        Receiver.__init__(self, Nx, Ny, dx, dy, Nt, dt)
        self.c = np.zeros([Nx, Ny])
        self.rho = np.zeros([Nx, Ny])

    def ReadVelocity(self, c):
        self.c = c

    def ReadDensity(self, rho):
        self.rho = rho

    def PrintModelInfo(self):
        print("------------------------------")
        print(f"Model grid size: Nx = {self.Nx}, Ny = {self.Ny}")
        print(f"Model spatial size: dx = {self.dx}, dy = {self.dy}")
        print(f"Time information: Nt = {self.Nt}, dt = {self.dt}")
        print(f"CFL number is: ")
        print("------------------------------")

class ModelPML(Model):

    def __init__(self, Nx, Ny, dx, dy, Nt, dt, pml_len, pml_alpha):
        Model.__init__(self, Nx, Ny, dx, dy, Nt, dt)
        self.pml_len = pml_len
        self.pml_alpha = pml_alpha
        self.Nx_pml = Nx + 2*pml_len
        self.Ny_pml = Ny + 2*pml_len
        self.c_pml = np.zeros([self.Nx_pml, self.Ny_pml])
        self.rho_pml = np.zeros([self.Nx_pml, self.Ny_pml])
        self.pml_value = []
        self.sigma_x = []
        self.sigma_y = []
        self.source_position_pml = []
        self.receiver_position_pml = []
        # self.SetPML()

    def SetPML(self, method=0):
        # method = 
        # 0: linear
        # 1: quadratic
        # 2: exponential, tbd

        if method == 0:
            self.pml_value = np.linspace(0, 1, self.pml_len)
            self.pml_value = self.pml_value * self.pml_alpha
        if method == 1:
            self.pml_value = np.linspace(0, 1, self.pml_len)
            self.pml_value = self.pml_value ** 2 * self.pml_alpha
        self.pml_value = np.flip(self.pml_value)

        self.Nx_pml = self.Nx + 2*self.pml_len
        self.Ny_pml = self.Ny + 2*self.pml_len
        self.c_pml = self.__ExpandPMLArea(self.c)
        self.rho_pml = self.__ExpandPMLArea(self.rho)
        self.source_position_pml = self.source_position + self.pml_len
        self.receiver_position_pml = self.receiver_position + self.pml_len
        # for i in range(self.source_num):
        #     self.source_position_pml.append(self.source_position[i] + self.pml_len)
        #     print(i)
        # for i in range(self.receiver_num):
        #     self.receiver_position_pml.append(self.receiver_position[i] + self.pml_len)
        #     print(i)

        self.sigma_x = np.zeros([self.Nx_pml, self.Ny_pml])
        self.sigma_y = np.zeros([self.Nx_pml, self.Ny_pml])
        for i in range(self.pml_len):
            self.sigma_x[i, :] = self.pml_value[i]
            self.sigma_x[-1-i, :] = self.pml_value[i]
            self.sigma_y[:, i] = self.pml_value[i]
            self.sigma_y[:, -1-i] = self.pml_value[i]
        self.sigma_x = self.sigma_x
        self.sigma_y = self.sigma_y

    def __ExpandPMLArea(self, A):
        A_pml = np.zeros([self.Nx_pml, self.Ny_pml])
        A_pml[self.pml_len : -self.pml_len, self.pml_len : -self.pml_len] = A
        for i in range(self.pml_len):
            A_pml[i, :] = A_pml[self.pml_len, :]
            A_pml[-1-i, :] = A_pml[-1-self.pml_len, :]
            A_pml[:,i] = A_pml[:,self.pml_len]
            A_pml[:,-1-i] = A_pml[:,-1-self.pml_len]
        return A_pml
    
    def PrintModelInfo(self):
        print("------------------------------")
        print(f"Model grid size: Nx = {self.Nx}, Ny = {self.Ny}")
        print(f"Model spatial size: dx = {self.dx}, dy = {self.dy}")
        print(f"Time information: Nt = {self.Nt}, dt = {self.dt}")
        print(f"CFL number is: ")
        print(f"PML information: pml_len = {self.pml_len}, pml_alpha = {self.pml_alpha}")
        print("------------------------------")

    def WriteModel(self):
        output_dict = {
            "Nx" : self.Nx,
            "Ny" : self.Ny,
            "dx" : self.dx,
            "dy" : self.dy,
            "Nt" : self.Nt,
            "dt" : self.dt,
            "pml_len" : self.pml_len,
            "pml_alpha" : self.pml_alpha,
            "c" : self.c.tolist(),
            "rho" : self.rho.tolist(),
            "source_num" : self.source_num,
            "source_position" : self.source_position.tolist(),
            "source_position_pml" : self.source_position_pml.tolist(),
            "source_fn" : self.source_fn.tolist(),
            "receiver_num" : self.receiver_num,
            "receiver_position" : self.receiver_position.tolist(),
            "receiver_position_pml" : self.receiver_position_pml.tolist(),
            "Nx_pml" : self.Nx_pml,
            "Ny_pml" : self.Ny_pml,
            "c_pml" : self.c_pml.tolist(),
            "rho_pml" : self.rho_pml.tolist(),
            "sigma_x" : self.sigma_x.tolist(),
            "sigma_y" : self.sigma_y.tolist()
        }

        print("Output to sample.json")
        # Serializing json
        json_object = json.dumps(output_dict, indent=4)
        with open("./data/temp_py.json", "w") as outfile:
            outfile.write(json_object)
        print("Output done.")
        print("------------------------------")


        
# def main():
#     Nx = 101
#     Ny = 101
#     dx = 0.01
#     dy = 0.01
#     Nt = 501
#     dt = 0.001

#     pml_len = 30
#     pml_alpha = 20

#     # c = np.random.randn(Nx, Ny)
#     # rho = np.random.randn(Nx, Ny)
#     # for i in range(Nx):
#     #     for j in range(Ny):
#     #         c[i,j] = i * Ny + j
#     #         rho[i,j] = i + j

#     c = np.ones([Nx, Ny])
#     rho = np.ones([Nx, Ny])

#     m = ModelPML(Nx, Ny, dx, dy, Nt, dt, pml_len, pml_alpha)

#     m.ReadVelocity(c)
#     m.ReadDensity(rho)
#     m.SetSource(3, [[i*10, 5] for i in range(3)], 5, 0.2)
#     m.SetReceiver(11, [[i*10, 0] for i in range(11)])
#     m.SetPML()

#     m.PrintModelInfo()

#     m.WriteModel()

#     # print(m.c)
#     # print(m.c_pml)
#     # print(m.sigma_x)
#     # print(m.source_position)
#     # print(m.source_position_pml)

#     print(m.source_fn.shape)




# if __name__ == "__main__":
#     main()
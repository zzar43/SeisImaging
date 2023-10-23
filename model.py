import numpy as np
import json

# define source function
class SourceData:
    def __init__(self):
        self.source_num = 0
        self.Nt = 0
        self.data = np.array([])

    def SetSource(self, source_num, Nt, dt, fre, center):
        self.source_num = source_num
        self.Nt = Nt
        t = np.linspace(0, (Nt-1)*dt, Nt)
        self.data = np.zeros([source_num, Nt])
        for i in range(source_num):
            self.data[i] = self.RickerFn(t, fre, center)

    def RickerFn(self, t, fre, center):
        return (1 - 2 * np.pi**2 * fre**2 * (t-center)**2) * np.exp(-1 * np.pi**2 * fre**2 * (t-center)**2)
    
    def WriteSource(self):
        json_dict = self.__dict__
        for key, value in self.__dict__.items():
            if type(value) == np.ndarray:
                json_dict[key] = value.flatten().tolist()

        with open('./data/source.json', 'w', encoding='utf-8') as f:
            json.dump(json_dict, f, ensure_ascii=False, indent=4)

class Mesh:
    def __init__(self, Nx, Ny, dx, dy):
        self.Nx = Nx
        self.Ny = Ny
        self.dx = dx
        self.dy = dy
        self.N = Nx * Ny

    def printMeshInfo(self):
        print("Model spatial size: Nx = %d, Ny = %d" % (self.Nx, self.Ny))
        print("Model spatial grid size: dx = %f, dy = %f" % (self.dx, self.dy))
        print("Model size: N = %d" % (self.N))

class Time:
    def __init__(self, Nt, dt):
        self.Nt = Nt
        self.dt = dt

    def printTimeInfo(self):
        print("Time info: Nt = %d, dt = %f" % (self.Nt, self.dt))

class Acquisition:
    def __init__(self, source_num, source_position, receiver_num, receiver_position, source_type=0):
        self.source_type = source_type
        self.source_num = source_num
        self.receiver_num = receiver_num

        if type(source_position) != np.ndarray:
            self.source_position = np.array(source_position)
        else:
            self.source_position = source_position

        if type(receiver_position) != np.ndarray:
            self.receiver_position = np.array(receiver_position)
        else:
            self.receiver_position = receiver_position

        if self.source_position.shape != (self.source_num, 2):
            raise IndexError("Please check dimension of source.")
        
        if self.receiver_position.shape != (self.receiver_num, 2):
            raise IndexError("Please check dimension of receiver.")

class Field:
    def __init__(self, M):
        self.data = np.zeros([M.Nx, M.Ny])

    def readData(self, data):
        if self.data.shape != data.shape:
            raise IndexError("Please check field shape and mesh shape.")
        self.data = data

# This model is for forward modelling and inversion
class Model:
    def __init__(self, m: Mesh, t: Time, c: np.ndarray, rho: np.ndarray, a: Acquisition):
        
        self.Nx, self.Ny = m.Nx, m.Ny
        self.dx, self.dy = m.dx, m.dy
        self.Nt, self.dt = t.Nt, t.dt

        self.c = c
        self.rho = rho

        self.source_num, self.source_position = a.source_num, a.source_position
        self.receiver_num, self.receiver_position = a.receiver_num, a.receiver_position

        # inverse problem
        self.c_inv = c
        self.rho_inv = rho

class ModelPML(Model):
    def __init__(self, m: Mesh, t: Time, c: np.ndarray, rho: np.ndarray, a: Acquisition, pml_len: int, pml_alpha: float, pml_method=0):

        Model.__init__(self, m, t, c, rho, a)

        self.pml_len = pml_len
        self.pml_alpha = pml_alpha
        
        self.Nx_pml = m.Nx + 2 * pml_len
        self.Ny_pml = m.Ny + 2 * pml_len

        self.source_position_pml = self.source_position + pml_len
        self.receiver_position_pml = self.receiver_position + pml_len

        self.c_pml = self.__ExpandPMLArea(self.c)
        self.rho_pml = self.__ExpandPMLArea(self.rho)

        self.__PMLMethod(pml_method)
        self.__PMLSigma()

    def __ExpandPMLArea(self, A):
        A_pml = np.zeros([self.Nx_pml, self.Ny_pml])
        A_pml[self.pml_len : -self.pml_len, self.pml_len : -self.pml_len] = A
        for i in range(self.pml_len):
            A_pml[i, :] = A_pml[self.pml_len, :]
            A_pml[-1-i, :] = A_pml[-1-self.pml_len, :]
            A_pml[:,i] = A_pml[:,self.pml_len]
            A_pml[:,-1-i] = A_pml[:,-1-self.pml_len]
        return A_pml
    
    def __PMLMethod(self, pml_method):
        if pml_method == 0:
            self.pml_value = np.linspace(0, 1, self.pml_len)
            self.pml_value = self.pml_value * self.pml_alpha
        if pml_method == 1:
            self.pml_value = np.linspace(0, 1, self.pml_len)
            self.pml_value = self.pml_value ** 2 * self.pml_alpha
        self.pml_value = np.flip(self.pml_value)

    def __PMLSigma(self):
        self.sigma_x = np.zeros([self.Nx_pml, self.Ny_pml])
        self.sigma_y = np.zeros([self.Nx_pml, self.Ny_pml])
        for i in range(self.pml_len):
            self.sigma_x[i, :] = self.pml_value[i]
            self.sigma_x[-1-i, :] = self.pml_value[i]
            self.sigma_y[:, i] = self.pml_value[i]
            self.sigma_y[:, -1-i] = self.pml_value[i]
        self.sigma_x = self.sigma_x
        self.sigma_y = self.sigma_y

    def WriteModel(self):
        json_dict = self.__dict__
        for key, value in self.__dict__.items():
            if type(value) == np.ndarray:
                json_dict[key] = value.flatten().tolist()

        with open('./data/data.json', 'w', encoding='utf-8') as f:
            json.dump(json_dict, f, ensure_ascii=False, indent=4)

# read seis data
def ReadWavefield(model):
    Nx = model.Nx
    Ny = model.Ny
    pml_len = model.pml_len

    f = open('./data/seis_data.json')
    data = json.load(f)
    f.close()
    
    last_wavefield = np.array(data["last_wavefield"]).reshape(Nx+2*pml_len, Ny+2*pml_len)
    return last_wavefield[pml_len:Nx+pml_len, pml_len:Ny+pml_len]

def ReadSeisMulti(model):
    source_num = model.source_num
    receiver_num = model.receiver_num
    Nt = model.Nt

    f = open('./data/seis_data.json')
    data = json.load(f)
    f.close()
    
    seis_signal_multi = np.array(data["seis_signal_multi"]).reshape(source_num, receiver_num, Nt)
    return seis_signal_multi

if __name__ == "__main__":

    Nx = 81
    Ny = 101
    Nt = 301
    dt = 0.001
    dx = 0.01
    dy = 0.01
    
    source_num = 9
    receiver_num = 101

    pml_len = 30
    pml_alpha = 20

    source_fre = 5
    source_center = 0.2

    m = Mesh(Nx, Ny, dx, dy)
    t = Time(Nt, dt)
    a = Acquisition(source_num, [[i*10,5] for i in range(source_num)],
                    receiver_num, [[i,0] for i in range(receiver_num)])
    
    c = np.ones([Nx, Ny])
    rho = np.ones([Nx, Ny])

    model = ModelPML(m, t, c, rho, a, pml_len, pml_alpha)
    model.WriteModel()

    s = SourceData()
    s.SetSource(source_num, Nt, dt, source_fre, source_center)
    s.WriteSource()
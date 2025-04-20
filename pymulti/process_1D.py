import numpy as np

class DataProcessor:
    def __init__(self):
        self.dataListAll = []
        self.dataArrayS2T = None
        self.dataArrayS1T = None

    def read_fort10_file(self, filename):
        """读取fort.10文件并提取dataArrayS2T数据"""
        with open(filename, 'r') as fid:
            self.dataListAll = [x.strip() for x in fid.readlines()]

        countDataAll = len(self.dataListAll)
        countS1 = int(self.dataListAll[0])
        countS2 = int(self.dataListAll[countS1 * 3 + 1])

        iBeginS1 = 1
        iEndS1 = countS1 * 3
        iBeginS2 = iEndS1 + 2
        iEndS2 = countDataAll

        columnsS1 = 3
        columnsS2 = (countDataAll - 1 - countS1 * 3 - 1) // countS2

        dataListS1 = self.dataListAll[iBeginS1:iEndS1]
        dataListS2 = self.dataListAll[iBeginS2:iEndS2]
        dataListS1.append(0)
        dataArrayS1 = np.reshape(dataListS1, (columnsS1, countS1))
        dataArrayS2 = np.reshape(dataListS2, (columnsS2, countS2))
        self.dataArrayS2T = dataArrayS2.T
        self.dataArrayS1T = dataArrayS1.T

    def extract(self, variable_name):
        """根据变量名提取dataArrayS2T中的相应变量，并进行平滑处理"""
        arrays_dict={}
        var_index=np.where(self.dataArrayS1T == 'CMC')[0]
        var_data = self.dataArrayS1T[var_index, 2:].astype(float)
        arrays_dict['M']=var_data
        for name in variable_name:
            if name == "V":
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                weights = np.array([0.5, 0.5])
                smoothed_data = np.zeros([1, var_data.shape[1]])
                for i in range(var_data.shape[1]):
                    rho = self.dataArrayS2T[np.where(self.dataArrayS2T == 'R')[0], 2:].astype(float)
                    RHO_smoothed = np.convolve(rho[:, i], weights, mode='valid')
                    RHO_smoothed = np.concatenate(([rho[0, i]], RHO_smoothed, [rho[-1, i]]))
                    var_data_negative = var_data[:, i] < 0
                    RHO_smoothed_negative = RHO_smoothed[var_data_negative]
                    if np.any(var_data_negative):
                        smoothed_data[0, i] = np.sum(RHO_smoothed_negative * var_data[var_data_negative, i]) / np.sum(
                            RHO_smoothed_negative)
                    else:
                        smoothed_data[0, i] = 0
                arrays_dict[name] = -smoothed_data/1e5,var_data/1e5

            elif name=="TI":
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                smoothed_data = np.zeros([1, var_data.shape[1]])
                M=arrays_dict['M'].reshape(-1)
                for i in range(var_data.shape[1]):
                    smoothed_data[0, i] = np.sum(M * var_data[:, i]) / np.sum(M)
                arrays_dict[name]=smoothed_data/1e3

            elif name=='VIMPLO':
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                arrays_dict[name]=var_data/1e5

            elif name=='IFAR':
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                arrays_dict[name]=var_data

            elif name=='RHORDT':
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                arrays_dict[name]=var_data

            elif name=="X":
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                arrays_dict[name]=var_data

            elif name=="TIME":
                var_index = np.where(self.dataArrayS2T == name)[0]
                if len(var_index) == 0:
                    raise ValueError(f"Variable '{name}' not found in dataArrayS2T.")
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                arrays_dict[name]=var_data

            elif name=='EKIMP':
                var_index = np.where(self.dataArrayS2T == 'V')[0]
                var_data = self.dataArrayS2T[var_index, 2:].astype(float)
                weights = np.array([0.5, 0.5])
                M_smoothed = np.convolve(arrays_dict['M'].reshape(-1), weights, mode='valid')
                M_smoothed = np.concatenate(([arrays_dict['M'][0,0]], M_smoothed, [arrays_dict['M'][-1,0]]))
                smoothed_data = np.zeros([1, var_data.shape[1]])
                for i in range(var_data.shape[1]):
                    var_data_negative = var_data[:, i] < 0
                    M_smoothed_negative = M_smoothed[var_data_negative]
                    if np.any(var_data_negative):
                        smoothed_data[0, i] = np.sum(M_smoothed_negative * var_data[var_data_negative, i]**2)/2
                    else:
                        smoothed_data[0, i] = 0
                arrays_dict['EKIMP'] = smoothed_data/5.6/1e10
        return arrays_dict


# 使用示例
# processor = DataProcessor()
# processor.read_fort10_file('fort_0.10')
# smoothed_neutrons = processor.extract('NEUTRONS')
# print(smoothed_neutrons)

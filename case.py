import indices

# 高度列表
height = [0, 10, 25, 50, 75, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490,
          520, 560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1120, 1160, 1200]

# 压强列表
pressure = [1008, 1005, 1000, 992, 974, 952, 925, 913, 878, 850, 843, 815, 783, 768, 700, 623, 606, 593, 582, 552, 549,
            532, 500, 448, 444, 441, 400, 374, 331, 323, 300, 287, 250]

# 温度列表
temperature = [15.6, 16.0, 15.8, 15.6, 14.4, 13.4, 11.6, 10.6, 7.6, 6.0, 5.6, 3.6, 2.2, 2.0, -3.9, -11.7, -13.1, -13.5,
               -12.1, -15.1, -14.9, -16.3, -20.3, -26.9, -26.5, -26.9, -32.7, -36.5, -43.9, -45.1, -49.5, -52.1, -57.9]

# 露点列表
dew_point = [-1.4, -2.0, -2.2, -3.4, -3.6, -6.6, -8.4, -5.4, -5.4, -4.0, -5.4, -6.4, -13.8, -17.0, -20.9, -20.7, -20.1,
             -21.5, -40.1, -33.1, -29.9, -32.3, -31.3, -35.9, -37.5, -39.9, -46.7, -50.5, -52.9, -53.1, -57.5, -59.1,
             -66.9]

# 风速列表（9999为缺失值）
w_speed = [2, 9999, 5, 9999, 9999, 9999, 7, 9999, 9999, 13, 9999, 9999, 9999, 9999, 17, 9999, 9999, 9999, 9999, 9999,
           9999, 9999, 34, 9999, 9999, 9999, 32, 9999, 9999, 9999, 45, 9999, 46]

# 风向列表
w_direct = [185, 9999, 190, 9999, 9999, 9999, 190, 9999, 9999, 230, 9999, 9999, 9999, 9999, 275, 9999, 9999, 9999, 9999,
            9999, 9999, 9999, 295, 9999, 9999, 9999, 310, 9999, 9999, 9999, 315, 9999, 305]

# 对流有效位能
CAPE = indices.CAPE_index(pressure, temperature, dew_point)
print('CAPE', CAPE)

# 全总指数
TT = indices.TT_Index(pressure, temperature, dew_point)
print('TT', TT)

# 杰弗森指数
JI = indices.JI_Index(pressure, temperature, dew_point)
print('JI', JI)

# A指数
A = indices.A_index(pressure, temperature, dew_point)
print('A', A)

# K指数
K = indices.K_index(pressure, temperature, dew_point)
print('K', K)

# S指数
S = indices.S_index(pressure, temperature, dew_point)
print('S', S)

# 抬升凝结高度
LCL = indices.LCL_index(pressure, temperature, dew_point)
print('LCL', LCL)

# 沙氏指数
SI = indices.SI_index(pressure, temperature, dew_point)
print('SI', SI)

# 条件性稳定度指数
IC = indices.IC_index(pressure, temperature, dew_point)
print('IC', IC)

# 抬升指数
LI = indices.LI_index(pressure, temperature, dew_point)
print('LI', LI)

# 深对流指数
DCI = indices.DCI_index(pressure, temperature, dew_point)
print('DCI', DCI)

# 自由对流高度
LFC = indices.LFC_index(pressure, temperature, dew_point)
print('LFC', LFC)

# 强天气威胁指数
SWEAT = indices.SWEAT_index(pressure, temperature, dew_point, w_speed, w_direct)
print('SWEAT', SWEAT)

# 逆温层
inv_hei = indices.inver_height(height, temperature)
print(inv_hei)

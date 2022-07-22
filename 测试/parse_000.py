import os.path

from indices import *
import pandas as pd

PATH = r"D:\文档\微波辐射计\北京指数\3月7-14日TlogP"

def get_stations(time):
    stations = []
    filename = str(time) + '.000'
    with open(os.path.join(PATH, filename), 'r', encoding='gbk') as f:
        for line in f:
            if line.startswith('5'):
                stations.append(line[:5])
    return stations

def get_data(station_id, time):
    """提取传入站点的数据，保存为temp文件"""
    station_id = str(station_id)
    time = str(time) + '.000'
    with open(os.path.join(PATH, time), 'r', encoding='gbk') as f:
        for line in f:
            if not line.startswith(' ') and line.startswith(station_id):
                with open('temp.txt', 'w') as temp:
                    temp.write(line)
                    temp.write('   pressure   height   temperature   dew_point   w_direct   w_speed\n')
                    line = f.readline()
                    while line.startswith(' '):
                        temp.write(line)
                        line = f.readline()
                    break

def to_df():
    data = pd.read_table('temp.txt', skiprows=1, sep=r'\s+')
    data.drop_duplicates(subset='pressure', keep='first', inplace=True)
    data = data[data["pressure"] > 200]
    return data


def cal_index(time, station):
    get_data(station, time)
    df = to_df()
    print('station', station, 'time', time)

    # CAPE_value = metpyTest.test_CAPE(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist()
    # )
    # print('CAPE', CAPE_value)

    # # 计算指数
    # CAPE_value = cape(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('CAPE', CAPE_value)

    # TT = TT_Index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('TT', TT)
    #
    # JI = JI_Index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('JI', JI)
    #
    # A = A_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('A', A)
    #
    # K = K_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('K', K)
    #
    # SI = SI_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('SI', SI)

    # LI = LI_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('LI', LI)

    # S = S_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('S', S)
    #
    # TQ = TQ_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('TQ', TQ)
    #
    # CT = CT_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('CT', CT)
    #
    # VT = VT_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    # )
    # print('VT', VT)
    #
    # KO = KO_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('KO', KO)
    #
    # LCL = LCL_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('LCL', LCL)
    #
    # IC = IC_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('IC', IC)
    #
    # BIC = BIC_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('BIC', BIC)
    #
    # DCI = DCI_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist(),
    # )
    # print('DCI', DCI)

    # LFC = LFC_index(
    #     df['pressure'].tolist(),
    #     df['temperature'].tolist(),
    #     df['dew_point'].tolist()
    # )
    # print("LFC", LFC)

    SWEAT = SWEAT_index(
        df['pressure'].tolist(),
        df['temperature'].tolist(),
        df['dew_point'].tolist(),
        df['w_speed'].tolist(),
        df['w_direct'].tolist(),
    )
    print("SWEAT", SWEAT)


if __name__ == '__main__':
    time = 22031420
    # time = 20210702200000
    stations = []
    filename = str(time) + '.000'
    # with open(os.path.join(PATH, filename), 'r', encoding='gbk') as f:
    #     for line in f:
    #         if line.startswith('5'):
    #             stations.append(line[:5])
    #
    # for station in stations:
    #     cal_index(
    #         time=time,
    #         station=station
    #     )
    cal_index(
        time=time,
        station=54511
    )
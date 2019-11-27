from pycwr.io.auto_io import radar_io
import pyart
import sys
import os

def save_cfradial(china_radar_file, save_file=None):
    """
    :param china_radar_file: radar data filename
    :param save_file: savename of cfradial format data
    :return:
    """
    radar = radar_io(china_radar_file).ToPyartRadar()
    if save_file is None:
        save_file = china_radar_file + ".nc"
    pyart.io.write_cfradial(save_file, radar)
    return 0

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("warning using!!! example: transfrom2cfradial filename  savename")
    elif not os.path.exists(sys.argv[1]):
        print("file is not exist!!!")
    elif len(sys.argv) == 2:
        save_cfradial(sys.argv[1])
    else:
        save_cfradial(sys.argv[1], sys.argv[2])

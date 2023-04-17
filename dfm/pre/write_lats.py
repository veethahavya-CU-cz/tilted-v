#!/usr/bin/python
import pickle


with open('../../common/data/ugrid.pkl', 'rb') as ugrid_file:
    ugrid_dict = pickle.load(ugrid_file)
ugrid_file.close()


ex = ugrid_dict['ex']
ey = ugrid_dict['ey']
with open('../lat.ext', 'w') as f:
    f.write("[General]\n")
    f.write("fileVersion = 2.01\n")
    f.write("fileType = extForce\n")
    for i, (x,y) in enumerate(zip(ex,ey)):
        f.write("\n[Lateral]\n")
        f.write("id = {:d}\n".format(i))
        f.write("name = lat_{:d}\n".format(i))
        f.write("numCoordinates = 1\n")
        f.write("xCoordinates = {:f}\n".format(x))
        f.write("yCoordinates = {:f}\n".format(y))
        f.write("discharge = realtime\n")
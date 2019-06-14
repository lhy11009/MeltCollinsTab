import os
import numpy as np
from ReadFiles import read_vtk_file
import Mathf as Mathf

_year = 365*24*3600


class NoFileError(Exception):
    pass


class FlagsError(Exception):
    pass

def bound_points(coord0, noz, nnoz = 1):
    coord1 = np.zeros((coord0.shape[0]//noz, coord0.shape[1]))
    for row in range(coord0.shape[0]):
        if row % noz == nnoz-1:
            nxy = row//noz
            coord1[nxy,:] = coord0[row,:]
    return coord1

def bound_elements(coord0,nox,noz):
    coord1 = np.zeros(((nox-1)**2,3))
    n0 = np.zeros(4)
    for i in range(nox-1):
        for j in range(nox-1):
            n0[0] = i*noz*nox+j*noz
            n0[1] = (i+1)*noz*nox+j*noz
            n0[2] = i*noz*nox+(j+1)*noz
            n0[3] = (i+1)*noz*nox+(j+1)*noz
            n = i*(nox-1)+j
            for k in n0:
                coord1[n,:] = coord1[n,:]+coord0[int(k),:]
            coord1[n,:] = coord1[n,:]/4.0
    return coord1

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


def adjustThetaPhi(coord, index):
    for i in range(coord.shape[0]):
        coord[i, index[0]] = np.pi/2.0 - coord[i, index[0]]
        if coord[i, index[1]] < 0.0:
            coord[i, index[1]] = coord[i, index[1]] + 2.0 * np.pi


def adjustThetaPhi1(coord, index):
    for i in range(coord.shape[0]):
        coord[i, index[0]] = np.pi/2.0 - coord[i, index[0]]
        if coord[i, index[1]] < -np.pi:
            coord[i, index[1]] += 2.0 * np.pi
        elif coord[i, index[1]] > np.pi:
            coord[i, index[1]] -= 2.0 * np.pi


def expandLonBound(_coord, _field, index, bd, limit):
    n = 0
    coord1 = np.zeros(_coord.shape)
    field1 = np.zeros(_field.shape)
    for i in range(_coord.shape[0]):
        lon = _coord[i, index]
        if lon < bd[0] + limit:
            coord1[n, :] = _coord[i, :]
            coord1[n, index] = _coord[i, index] + 2.0 * np.pi
            field1[n] = _field[i]
            n = n + 1
        elif lon > bd[1] - limit:
            coord1[n, :] = _coord[i, :]
            coord1[n, index] = _coord[i, index] - 2.0 * np.pi
            field1[n] = _field[i]
            n = n + 1
    coord = np.concatenate((_coord, coord1[0:n, :]), axis=0)
    field = np.concatenate((_field, field1[0:n, :]), axis=0)
    return coord, field


class CombineSurfFile():
    def __init__(self, fname, _inD, **kwargs):
        self.fname = fname
        nprocx = int(_inD['nprocx'])
        nprocz = int(_inD['nprocz'])
        self.nprocz = nprocz
        nodex = int(_inD['nodex'])
        nodez = int(_inD['nodez'])
        self.nproc = 12 * nprocx * nprocx * nprocz
        self.nox = (nodex+1)//nprocx
        self.noz = (nodez+1)//nprocz
        try:
            nodez = kwargs['nodez']
        except KeyError:
            nodez = 1
        if nodez < 2:
            self.procz = 0
            self.nnoz = nodez
        else:
            self.procz = (nodez-2)//(self.noz-1)
            self.nnoz = 2 + (nodez-2) % (self.noz-1)

    def __call__(self, route, case_name, step, procz, flags, **kwargs):
        print("Combine surf file, step = %s" % step)
        try:
            _oroute = kwargs['oroute']
        except KeyError:
            _oroute = route
        try:
            _area = kwargs['area']
        except KeyError:
            _area = False
        try:
            _expandLon = kwargs['expandLon']
        except KeyError:
            _expandLon = False
        if flags == "layer":
            self.CallLayer(route, case_name, step, procz)
        elif flags == "surf":
            self.CallSurf(route, case_name, step, procz,
                          obl=kwargs["obl"], expandLon=_expandLon, oroute=_oroute,
                          area=_area)
        else:
            raise FlagsError("flags %s is wrong" % (flags))

    def CallLayer(self, route, case_name, step, procz):
        nprocxy = self.nproc//self.nprocz
        out_filename = os.path.join(route, "%s.lgmt.%s.%d.%d" %
                                           (case_name, self.fname, self.nodez, step))
        mm0 = np.zeros(2)
        print(out_filename)
        if os.path.isfile(out_filename):
            os.remove(out_filename)
        for i in range(nprocxy):
            field0, coord0 = self.Read1(route, case_name, step, i, self.procz)
            coord1 = bound_points(coord0, self.noz)
            field1 = bound_points(field0, self.noz, self.nnoz)
            mm1 = self.MaxMin(field1, 0)
            mm0 = mm1
            coord1 = Mathf.cart2sph(coord1)
            adjustThetaPhi(coord1, [1, 2])
            self.Output(out_filename, field1, coord1)
        line_prepender(out_filename, "%s %s" % (mm0[0], mm0[1]))
        pass

    def CallSurf(self, route, case_name, step, procz, **kwargs):
        try:
            oroute = kwargs['oroute']
        except KeyError:
            oroute = route
        nprocxy = self.nproc//self.nprocz
        out_filename = os.path.join(oroute, "%s.MF0.%d" % (case_name, step))
        try:
            expandLon = kwargs['expandLon']
        except KeyError:
            expandLon = False
        if kwargs['obl'] is False and os.path.isfile(out_filename):
            print("%s exists" % out_filename)
            return
        mm0 = np.zeros(2)
        if os.path.isfile(out_filename):
            os.remove(out_filename)
        #   append integration as header
        try:
            _area = kwargs['area']
        except KeyError:
            _area = False
        if _area is not False:
            _getInteg = True
            _integrate = 0.0
        else:
            _getInteg = False
        for i in range(nprocxy):
            MF0, coord0 = self.Read(route, case_name, step, i)
            coord1 = bound_elements(coord0, self.nox, self.noz)
            mm1 = self.MaxMin(MF0, 0)
            mm0[0] = min(mm0[0], mm1[0])
            mm0[1] = max(mm0[1], mm1[1])
            coord1 = Mathf.cart2sph(coord1)
            adjustThetaPhi1(coord1, [1, 2])
            if expandLon:
                coord1, MF0 = expandLonBound(coord1, MF0, 2, (-np.pi, np.pi), np.pi/40.0)
            if i == 0:
                self.Output(out_filename, MF0, coord1)
            else:
                self.Output(out_filename, MF0, coord1, append=True)
            if _getInteg:
                _integrate1 = self.IntegrateField(MF0, _area)
                _integrate = _integrate + _integrate1
        _line = "%s %s" % (mm0[0], mm0[1])
        if _getInteg:
            for i in range(_integrate.size):
                _line += "  %.4e" % _integrate[i]
        line_prepender(out_filename, _line)
        pass

    #   .MF_Total file
    def TotalAmount(self, _cls, stepTuple, timeArray, **kwargs):
        print("TotalAmount")
        try:
            obl = kwargs['obl']
        except KeyError:
            obl = False
        route = _cls.route
        caseName = _cls.caseName
        timeScale = _cls.timeScale
        tScaleYear = timeScale/_year
        oFile = os.path.join(route, "%s.MF_Total" % caseName)
        if obl is False and os.path.isfile(oFile):
            return
        inFile = os.path.join(route, "%s.MF0.%d" % (caseName, stepTuple[0]))
        if os.path.isfile(inFile) is False:
            raise NoFileError(inFile)
        # Read in coordinate
        inData = np.genfromtxt(inFile, skip_header=1)
        shape = inData.shape
        coord0 = np.zeros((shape[0], 3))
        MM = np.zeros((shape[0], shape[1] - 2))
        coord0[:, 1:3] = inData[:, 0: 2]
        tt0 = timeArray[0]
        # melt rate * time
        for step in stepTuple:
            print("MF Total step: %d" % step)
            tt1 = timeArray[step]
            tt = tt1 - tt0
            tt0 = tt1
            inFile = os.path.join(route, "%s.MF0.%d" % (caseName, step))
            inData = np.genfromtxt(inFile, skip_header=1)
            MF = inData[:, 2: shape[1] + 1]
            MM = MM + MF * tt * tScaleYear
        mm0 = self.MaxMin(MM, 0)
        self.Output(oFile, MM, coord0)
        _line = "%s %s" % (mm0[0], mm0[1])
        try:
            _area = kwargs['area']
        except KeyError:
            pass
        else:
            _integrate = self.IntegrateField(MM, _area)
            for i in range(_integrate.size):
                _line += " %.4e" % (_integrate[i])
        line_prepender(oFile, _line)

    #   Intergration by element
    def IntegrateField(self, _field, _area):
        integrate = np.zeros(_field.shape[1])
        for i in range(_field.shape[0]):
            integrate = integrate + _field[i, :] * _area
        return integrate

    # Read in MF file, flag is unused
    def Read(self, route, case_name, step, procxy, flag="add"):
        proc = procxy*self.nprocz
        filename=os.path.join(route,"%s.MF.%d.%d"%(case_name,proc,step))
        MF0 = np.genfromtxt(filename)
        for j in range(1,self.nprocz):
            proc = procxy*self.nprocz+j
            filename = os.path.join(route,"%s.proc%d.0.vts"%(case_name,proc))
            ddict,coord = read_vtk_file(filename)
            filename=os.path.join(route,"%s.MF.%d.%d"%(case_name,proc,step))
            MF1 = np.genfromtxt(filename)
            MF0 = MF0 + MF1
        return MF0,coord

    def Read1(self,route,case_name,step,procxy,procz):
    #Read in MF file, flag is unused#
        proc = procxy*self.nprocz+procz
        filename=os.path.join(route,"%s.proc%d.%d.vts"%(case_name,proc,step))
        ddict,coord = read_vtk_file(filename)
        field = ddict[self.fname]
        return field,coord

    def OutputHeader(selt,filename,mm):
        with open(filename,'w') as f:
            f.write("%.4e %.4e\n"%(mm[0],mm[1]))

    # output to a MF0 file for gmt
    def Output(self, filename, MF0, coord, **kwargs):
        cshape = coord.shape
        Mshape = MF0.shape
        try:
            _append = kwargs['append']
        except KeyError:
            _append = False
        if _append is True:
            _type = 'a'
        else:
            _type = 'w'
        with open(filename, _type) as f:
            for i in range(Mshape[0]):
                f.write("%.5e" % coord[i, 1])
                for j in range(2, cshape[1]):
                    f.write(" %.5e" % coord[i, j])
                for j in range(Mshape[1]):
                    f.write(" %.5e" % MF0[i, j])
                f.write("\n")

    # check data range
    def CheckRange(self, filename, MF0, coord, **kwargs):
        pass

    def MaxMin(self,MF,column):
        mm = np.zeros(2)
        mm[0] = np.min(MF[:,column])
        mm[1] = np.max(MF[:,column])
        return mm


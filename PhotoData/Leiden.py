import numpy as np
import os

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

def Leiden(sp):
    fil = open(rootdir+'Leiden/all_cross_sections_text_continuum_0.1nm/'+sp+'_0.1nm.txt')
    lines = fil.readlines()
    fil.close()

    ind = np.nan
    for i,line in enumerate(lines):
        if "# wavelength   photoabsorption" in line:
            ind = i
            break

    keys = lines[ind].split()[1:]
    tmp = {}
    for key in keys:
        tmp[key] = []
    for line in lines[ind+1:]:
        dat = [float(a) for a in line.split()]
        for i,key in enumerate(keys):
            tmp[key].append(dat[i])

    xs = np.array(tmp['photoabsorption'])
    qy = []
    for i in range(len(xs)):
        if xs[i] == 0:
            qy.append(0)
        else:
            qy.append(tmp['photodissociation'][i]/xs[i])
    qy = np.array(qy)

    fil = open(rootdir+'Leiden/cross_section_properties.csv','r')
    lines = fil.readlines()
    fil.close()
    tmp1 = {}
    for i,line in enumerate(lines[19:]):
        species = line.split(',')[0].strip()
        prod = line.split(',')[-1].strip()
        tmp1[species] = prod

    data = {}
    data['wavelength'] = np.array(tmp['wavelength'])
    data['cross section'] = xs
    data[tmp1[sp]] = qy
    return data
    
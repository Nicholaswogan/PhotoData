import re
import sys
import numpy as np
import os

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

class phidrates():
    
    def __init__(self,spec):
        fil = open(rootdir+'phidrates/phidrates.txt','r')
        lines = fil.readlines()
        fil.close()
        species = [[ln.strip() for ln in line.strip().split('|')][0] for line in lines]

        # find species
        try:
            spec_ind = species.index(spec)
        except:
            sys.exit(spec+' not in phidrates database')

        # get reactions
        if len(lines[spec_ind].split('|')) == 2:
            reactions = []
        else:
            reactions = [ln.strip() for ln in lines[spec_ind].split('|')][2:]

        # open data file
        fil = open(rootdir+'phidrates/'+spec+'.txt','r')
        file = fil.readlines()
        fil.close
        self.file = file

        inds = []
        for i,dat in enumerate(file):
            if len(dat)>0:
                if dat[0]=='0':
                    inds.append(i)
        self.inds = inds

        # get number of branches
        if not file[inds[-1]].split()[-1]=='branches':
            sys.exit()
        num_branches = int(file[inds[-1]].split()[-2])

        # branch names
        name_branches = file[inds[-1]+1].split()[1:]

        # check
        if len(name_branches) != num_branches+1 or len(name_branches) != len(reactions):
            sys.exit('here')

        # separate into ion and non-ion reactions
        neutral_inds = [i for i,name in enumerate(name_branches) if name.find('/')>-1 and name.find('+')==-1]
        ion_inds = [i for i in range(len(name_branches)) if i not in neutral_inds and name_branches[i] != 'Total']

        # get data
        data = {}
        if 'Total' in name_branches or 'total' in name_branches:
            # nm
            data['wavelength'] = np.array([float(dat.strip().split()[0])/10. for dat in file[inds[-1]+2:]])
            data['Total'] = np.array([float(dat.strip().split()[1]) for dat in file[inds[-1]+2:]])
            factors = np.array([])
            for i,dat in enumerate(file[inds[-1]+2:]):
                temp = np.sum([float(da) for da in dat.strip().split()[2:]])
                factor = data['Total'][i]/temp
                if factor > 1.1 or factor < 0.9:
                    print('Warning! Quantum yields not summing to 1 for '+spec)
                factors = np.append(factors,factor)
            for i,name in enumerate(name_branches[1:]):
                data[name] = np.array([float(dat.strip().split()[i+2]) for dat in file[inds[-1]+2:]])/data['Total']
                if np.max(data[name])>1.1:
                    print('Warning! Quantum yields not summing to 1 for '+spec)
        else:
            data['wavelength'] = np.array([float(dat.strip().split()[0])/10. for dat in file[inds[-1]+2:]])
            data['Total'] = np.array([float(dat.strip().split()[1]) for dat in file[inds[-1]+2:]])
            data[name_branches[0]] = np.array([float(dat.strip().split()[1]) for dat in file[inds[-1]+2:]])/data['Total']
            if len(name_branches) != 1:
                sys.exit('Issue parsing data')
        
        self.neutral = rx()
        self.ion = rx()
        
        # get neutral data
        self.neutral.branches = [name_branches[ind] for ind in neutral_inds]
        self.neutral.reactions = [reactions[ind] for ind in neutral_inds]
        self.neutral.data = {}
        self.neutral.data['wavelength'] = data['wavelength']
        self.neutral.data['cross section'] = data['Total']
        for name in self.neutral.branches:
            self.neutral.data[name] = data[name]

        # get ion data
        self.ion.branches = [name_branches[ind] for ind in ion_inds]
        self.ion.reactions = [reactions[ind] for ind in ion_inds]
        self.ion.data = {}
        self.ion.data['wavelength'] = data['wavelength']
        self.ion.data['cross section'] = data['Total']
        for name in self.ion.branches:
            self.ion.data[name] = data[name]
        self.get_meta_data()
            
            
    def get_meta_data(self):
        file = self.file
        inds = self.inds
        meta_data = {}
        for j in range(len(inds)-1):
            ref_temp = file[inds[j]:inds[j+1]]
            refs_for = ref_temp[0].strip().split()[-1]
            ref_ind = []
            for i in range(1,len(ref_temp)):
                if len(ref_temp[i].strip().split()) != 0:
                    if ref_temp[i].strip().split()[0] == 'Lambda':
                        ref_ind.append(i)
            ref_ind.append(len(ref_temp)-1)

            # get refs
            wv_range = []
            rough_ref = []
            for i in range(len(ref_ind)-1):
                try:
                    temp = ref_temp[ref_ind[i]].strip().split()
                    if temp[3] != 'A':
                        sys.exit('Issue parsing refs for '+spec)
                    # nm
                    if ref_temp[ref_ind[i]+1].strip().split()[0]=='band':
                        sys.exit()
                    wv_range.append([float(re.sub("[^0123456789\.]",'',a))/10. for a in temp[2].replace('-',' ').split()])
                    # get references
                    temp1 = [ref for ref in ref_temp[ref_ind[i]+1:ref_ind[i+1]]]
                    temp2 = []
                    for tt in temp1:
                        if tt[:3] != '   ':
                            break
                        if tt.strip()[0]=='(':
                            pass
                        elif tt.strip().split()[0]=='band':
                            pass
                        else:
                            temp2.append(tt.strip().replace('?',''))
                    each_ref = []
                    for i,tt in enumerate(temp2):
                        try:
                            float(tt[tt.find('(')+1:tt.find(')')])
                            each_ref.append(i)
                        except:
                            try:
                                [float(a) for a in tt[tt.find('(')+1:tt.find(')')].split(',')]
                                each_ref.append(i)
                            except:
                                pass
                    each_ref.append(len(temp2))
                    if len(each_ref)>1:
                        temp3 = [' '.join(temp2[each_ref[i]:each_ref[i+1]]) for i in range(len(each_ref)-1)]   
                        rough_ref.append(temp3)
                    else:
                        rough_ref.append(temp2)
                except:
                    rough_ref = rough_ref[:len(wv_range)]
                    pass
            temp_meta = {}
            temp_meta['wavelengths'] = wv_range
            temp_meta['references'] = rough_ref
            meta_data[refs_for] = temp_meta
        
        # get bibtex citations for refs
        self.loadbibs()
        for key in meta_data.keys():
            meta = meta_data[key]['references']
            temp1 = []
            for met in meta:
                temp2 = []
                for me in met:
                    temp2.append(self.get_bib(me))
                temp1.append(temp2)
            meta_data[key]['bibtex'] = temp1
        self.meta_data = meta_data
    
    def loadbibs(self):
        fil = open(rootdir+'phidrates/phidrates_refs.txt','r')
        lines_references = fil.readlines()
        fil.close()

        fil = open(rootdir+'phidrates/phidrates_bibs.txt','r')
        file_bibs = fil.readlines()
        fil.close()

        # re-organize bibs
        inds = []
        for i,line in enumerate(file_bibs):
            if line[0] == '@' or line=='No Google Scholar results.\n':
                inds.append(i)
        inds.append(len(file_bibs)-1)
        lines_bibs = []
        for i in range(len(inds)-1):
            lines_bibs.append("".join(file_bibs[inds[i]:inds[i+1]]))
        lines_bibs1 = []
        for i,line in enumerate(lines_bibs):
            if re.search(lines_references[i].strip().split(',')[0], line.split('\n')[0], re.IGNORECASE):
                lines_bibs1.append(line)
            else:
                lines_bibs1.append('No Google Scholar results.\n')
        # create dict
        self.ref2bib = {}
        for i,line in enumerate(lines_references[:-1]):
            self.ref2bib[line.strip('\n')] = lines_bibs[i]
    
    def get_bib(self,ref):
        try:
            bib = self.ref2bib[ref]
            if bib == 'No Google Scholar results.\n':
                sys.exit()
        except:
            bib = 'No bibtex citation avaliable'
        return bib
    
    
    def get_atmos_data(self,species):
        folder = rootdir+'XSECTIONS_alinc/'
        wv_vpl = []
        xs_vpl = []
        try:
            fil = open(folder+species+'/'+species+'.XS.dat','r')
            lines = fil.readlines()
            fil.close()
        except:
            sys.exit(species+' is not in the atmos database')
        try:
            if len(lines[2].split()) == 1 or len(lines[2].split()) == 2:
                jjj = 1
            else:
                jjj = lines[2].split().index('300K')

            for line in lines[4:]:
                wv_vpl.append(.1*float(line.split()[0]))
                xs_vpl.append(float(line.split()[jjj]))

            qy_files = [file for file in os.listdir(folder+species+'/') if file.endswith(".QY.dat")]
            prods = ['/'.join(qy.strip('.QY.dat').split('_')[2:]) for qy in qy_files]
            wvs = []
            qys = []
            for file in qy_files:
                fil = open(folder+species+'/'+file,'r')
                lines = fil.readlines()
                fil.close()
                wv = np.array([float(line.split()[0]) for line in lines[4:]])/10.
                qy = np.array([float(line.split()[1]) for line in lines[4:]])
                wvs.append(wv)
                qys.append(qy)

            start = np.min([wv.min() for wv in wvs])
            end = np.max([wv.max() for wv in wvs])
            resolution = .1
            wavs = np.arange(start,end,0.1)

            atmos_data = {}
            atmos_data['wavelength'] = wavs
            atmos_data['cross section'] = np.interp(wavs,wv_vpl,xs_vpl)
            for i,prod in enumerate(prods):
                atmos_data[prod] = np.interp(wavs,wvs[i],qys[i])
            self.atmos_data = atmos_data
            self.atmos_branches = prods
        except:
            print('Failed to get atmos data')
            self.atmos_data = {}
            pass
        
        
class rx():
    def __init__(self):
        self.reactions = []
        self.branches = []
        self.data = {}
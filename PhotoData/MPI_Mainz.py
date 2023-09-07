from bs4 import BeautifulSoup
import requests
import numpy as np
import sys
import re
import os
from doi2bib.crossref import get_bib
from matplotlib import pyplot as plt
from tabulate import tabulate
import subprocess

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'


class MPI_Mainz():

    def __init__(self,species):
        # options
        self.T_low = 0
        self.T_high = 600
        self.max_studies = np.inf
        self.bin_width = 10
        self.species = species
        self.verbose = True
        self.folder = rootdir+'XSECTIONS_alinc/'

        # search for species
        search = 'https://uv-vis-spectral-atlas-mainz.org/uvvis/search_species.html?csrfmiddlewaretoken=k6ilihVEfSjQXFJuVwD2vMmlihTQOkjR3scabAhJwxqs2XtQ1IZFQzkO0o1BUXMG&x=0&y=0&search='+species
        r1 = requests.get(search)
        html1 = r1.text
        ind1 = html1.find('results for')
        ind1 = html1[ind1:].find('<table')+ind1
        ind2 = html1[ind1:].find('</table>')+ind1
        searchp = BeautifulSoup(html1[ind1:ind2],"html.parser")
        jj = 0
        i = 0
        while i<1000:
            i+=1
            try:
                spec = searchp.find_all('tr')[i].a.text.strip()
                if spec == species:
                    jj = 1
                    self.spurl = searchp.find_all('tr')[i].a['href']
                    break
            except Exception:
                pass
        if jj==0:
            sys.exit('Species not found in MPI-Mainz database')

    def get_data(self,output = False):
        url = 'http://satellite.mpic.de/spectral_atlas/'+self.spurl
        r = requests.get(url)
        html = r.text
        ind1 = html.find('Data Sets:')
        ind2 = html[ind1:].find('</table>')+ind1
        parsed = BeautifulSoup(html[ind1:ind2],"html.parser")

        i = 0
        files = []
        papers = []
        yr = []
        temp = []
        for i in range(1,len(parsed.table.find_all('tr'))):
            details = parsed.table.find_all('tr')[i].find_all('td')

            # only consider data that is a function of wavelength
            if details[3].text.find('-')<0 and \
               details[-1].text.find('-')>0:
                try:
                    tempp = float(details[3].text[:-1])
                except ValueError:
                    tempp = 300
                if tempp>self.T_low and tempp<self.T_high:
                    try:
                        paper = parsed.table.find_all('tr')[i].find_all('td')[2].text
                        yr.append(float(paper[paper.find('(')+1:paper.find(')')]))
                        files.append(parsed.table.find_all('tr')[i].td.a['href'])
                        papers.append(paper)
                        temp.append(tempp)
                    except:
                        yr = yr[:len(temp)]
                        files = files[:len(temp)]
                        papers = papers[:len(temp)]
                        pass

        pre = 'http://joseba.mpch-mainz.mpg.de/spectral_atlas_data/'

        yr1 = []
        files1 = []
        papers1 = []
        refs = []
        temp1 = []
        wv = []
        xs = []
        for i in range(len(files)):
            if self.verbose:
                print(i,' files parsed',end='\r')
            r = requests.get(pre+files[i])
            data1 = r.text.replace('\t',' ').replace('\r', '').split('\n')
            wvv1 = []
            xss1 = []
            for dat in data1:
                try:
                    wvv, xss = [float(a) for a in dat.split()]
                    wvv1.append(wvv)
                    xss1.append(xss)
                except:
                    pass
            if len(wvv1)==0 or np.min(xss1)<=0:
                pass
            else:
                wv.append(wvv1)
                xs.append(xss1)
                yr1.append(yr[i])
                files1.append(files[i])
                papers1.append(papers[i])
                temp1.append(temp[i])

                # get a bibtex citation
                try:
                    pre_details = 'http://satellite.mpic.de/spectral_atlas/'
                    rr = requests.get(pre_details+files[i])
                    ind1 = [m.start() for m in re.finditer('DOI', rr.text)][1]
                    ind2 = rr.text[ind1:].find('</a>')+ind1+4
                    soup = BeautifulSoup(rr.text[ind1:ind2],"html.parser")
                    out = get_bib(soup.a['href'])
                    refs.append(out[1])
                except:
                    refs.append('Resource not found.')

                if len(wv)>self.max_studies:
                    break

        # build dict
        out = {}
        out['wavelength'] = wv
        out['cross section'] = xs
        out['temperature'] = temp1
        out['year'] = yr1
        out['papers'] = papers1
        out['bibtex'] = refs

        self.all_data = out
        
        # flatten data
        self.all_wv = np.array([item for w in wv for item in w])
        self.all_xs = np.array([item for w in xs for item in w])
        temps = np.array([])
        for i in range(len(temp1)):
            temps = np.append(temps,np.ones(len(wv[i]))*temp1[i])
        self.all_temps = temps

        if output:
            return out


    def find_best_data(self,wv_range = None, best = 'max resolution', output = False):

        wv = self.all_data['wavelength']
        xs = self.all_data['cross section']
        yr = self.all_data['year']
        papers = self.all_data['papers']
        refs = self.all_data['bibtex']

        rng = []
        for w in wv:
            rng.append([np.min(w),np.max(w)])
        tot_rng = [np.min(np.array(rng).flatten()),np.max(np.array(rng).flatten())]

        wavelengths = np.array([])
        cross_section = np.array([])
        citation = []
        bib = []

        if best == 'single longest':
            ind = np.argmax([rn[1]-rn[0] for rn in rng])

            wavelengths = np.append(wavelengths,np.array(wv[ind])[:])
            cross_section = np.append(cross_section,np.array(xs[ind])[:])
            citation.append(papers[ind])
            bib.append(refs[ind])
            citations = [citation[0],[wavelengths.min(),wavelengths.max()],bib[0]]
            citations1 = {}
            citations1['citation'] = [citations[0]]
            citations1['wavelength range'] = [citations[1]]
            citations1['bibtex'] = [citations[2]]

            best_data = {}
            best_data['wavelength'] = wavelengths
            best_data['cross section'] = cross_section

            self.best_data_citations = citations1
            self.best_data = best_data

            if output:
                return best_data, citations1
        else:

            if wv_range==None:
                bins = np.arange(tot_rng[0],tot_rng[1],self.bin_width)
                bins[-1] = tot_rng[1]
            else:
                bins = np.arange(wv_range[0],wv_range[1],self.bin_width)
                bins[-1] = wv_range[1]

            for j in range(len(bins)-1):
                options = []
                for i in range(len(rng)):
                    if rng[i][0]<=bins[j] and rng[i][1]>=bins[j+1] and np.min(xs[i])>=0:
                        temp_ind = np.where(((np.array(wv[i])<=bins[j+1]) & (bins[j]<=np.array(wv[i]))))[0]
                        if len(temp_ind)==0:
                            temp_res = 0
                        else:
                            temp_wv = np.array(wv[i])[temp_ind]
                            temp_xs = np.array(xs[i])[temp_ind]
                            if (np.max(temp_wv)-np.min(temp_wv))==0:
                                temp_res = 0
                            else:
                                temp_res = len(temp_wv)/(np.max(temp_wv)-np.min(temp_wv))
                        options.append([i,temp_res,yr[i],np.mean(np.log10(temp_xs))])
                if len(options)==0:
                    citation.append('No data')
                    bib.append('No data')
                    continue

                if best == 'max resolution':
                    max_res = np.argmax([val[1] for val in options])
                    ind = options[max_res][0]

                elif best == 'most recent':
                    max_yr = np.argmax([val[2] for val in options])
                    ind = options[max_yr][0]

                elif best == 'smallest xsections':
                    smallest = np.argmin([val[3] for val in options])
                    ind = options[smallest][0]

                elif best == 'largest xsections':
                    smallest = np.argmax([val[3] for val in options])
                    ind = options[smallest][0]

                inds = np.where(((np.array(wv[ind])<=bins[j+1]) & (bins[j]<=np.array(wv[ind]))))[0]
                wavelengths = np.append(wavelengths,np.array(wv[ind])[inds])
                cross_section = np.append(cross_section,np.array(xs[ind])[inds])
                citation.append(papers[ind])
                bib.append(refs[ind])

            citations = []
            j = 0
            bk_flag = False
            while True:
                for i in range(j,len(citation)):
                    if not citation[j]==citation[i]:
                        citations.append([citation[j],[bins[j],bins[i]],bib[j]])
                        j = i
                        break
                    if i == len(citation)-1:
                        citations.append([citation[j],[bins[j],bins[-1]],bib[j]])
                        bk_flag= True
                if bk_flag:
                    break
                if i>10000:
                    sys.exit('Something went wrong while sorting references')

            citations1 = {}
            citations1['citation'] = [cit[0] for cit in citations]
            citations1['wavelength range'] = [cit[1] for cit in citations]
            citations1['bibtex'] = [cit[2] for cit in citations]

            best_data = {}
            best_data['wavelength'] = wavelengths
            best_data['cross section'] = cross_section

            self.best_data_citations = citations1
            self.best_data = best_data

            if output:
                return best_data, citations1

    def get_atmos_data(self,species,output = False):
        wv_vpl = []
        xs_vpl = []
        try:
            fil = open(self.folder+species+'/'+species+'.XS.dat','r')
            lines = fil.readlines()
            fil.close()
            for line in lines[4:]:
                wv_vpl.append(.1*float(line.split()[0]))
                xs_vpl.append([float(a) for a in line.split()[1:]])
        except:
            pass

        wv_alin = []
        xs_alin = []
        try:
            file = self.folder+species+'/'+species+'_alinc.dat'
            fil = open(file,'r')
            lines = fil.readlines()
            fil.close()
            for line in lines[8:]:
                wv_alin.append(.1*float(line.split()[0]))
                xs_alin.append(float(line.split()[1]))
        except:
            pass

        wv_mpi = []
        xs_mpi = []
        try:
            file = self.folder+species+'/'+species+'_mpi.abs'
            fil = open(file,'r')
            lines = fil.readlines()
            fil.close()
            for line in lines[4:]:
                wv_mpi.append(.1*float(line.split()[0]))
                xs_mpi.append(float(line.split()[1]))
        except:
            pass

        wv_zahnle = []
        xs_zahnle = []

        try:
            file = self.folder+species+'/'+species+'_zahnle.abs'
            fil = open(file,'r')
            lines = fil.readlines()
            fil.close()
            for line in lines[2:]:
                wv_zahnle.append(.1*float(line.split()[0]))
                xs_zahnle.append(float(line.split()[1]))
        except:
            pass

        out = {}
        out['alinc'] = [wv_alin,xs_alin]
        out['vpl'] = [wv_vpl,xs_vpl]
        out['mpi'] = [wv_mpi,xs_mpi]
        out['zahnle'] = [wv_zahnle,xs_zahnle]

        self.atmos_data = out

        if output:
            return out


    def plot(self, plot_atmos = False, save = None):

        wv = self.all_data['wavelength']
        xs = self.all_data['cross section']
        temp = self.all_data['temperature']

        wv1 = self.all_wv
        xs1 = self.all_xs
        temps = self.all_temps

        plt.rcParams.update({'font.size': 15})
        if plot_atmos:
            fig,[ax,ax1] = plt.subplots(1,2,figsize=[17,4])
        else:
            fig,ax = plt.subplots(1,1,figsize=[7,4])

        map1 = ax.scatter(wv1,xs1,alpha=.5,c=temps,cmap='viridis',s=10)
        cbar = fig.colorbar(map1, ax=ax)
        cbar.set_label('Temperature (K)')
        ax.plot(self.best_data['wavelength'],self.best_data['cross section'],'r',label='Data Used')
        ax.set_yscale('log')
        ax.set_ylabel('Cross section (cm$^2$ molecule$^{-1}$)')
        ax.set_xlabel('Wavelength (nm)')
        ax.grid()
        # ax.set_xlim(100,400)
        # ax.set_ylim(1e-27,1e-15)
        # ax.legend()

        if plot_atmos:
            curr = self.atmos_data
            if len(curr['vpl'][0])>0:
                ax1.plot(curr['vpl'][0],curr['vpl'][1],'C0--',label=self.species+'.XS.dat')
            if len(curr['alinc'][0])>0:
                ax1.plot(curr['alinc'][0],curr['alinc'][1],'C1--',label=self.species+'_alinc.dat')
            if len(curr['mpi'][0])>0:
                ax1.plot(curr['mpi'][0],curr['mpi'][1],'C2--',label=self.species+'_mpi.abs')
            if len(curr['zahnle'][0])>0:
                ax1.plot(curr['zahnle'][0],curr['zahnle'][1],'C2--',label=self.species+'_zahnle.abs')
            ax1.set_yscale('log')
            ax1.set_xlabel('Wavelength (nm)')
            map1 = ax.scatter(1,-1,alpha=.5,c=[1],cmap='viridis',s=10)
            cbar = fig.colorbar(map1, ax=ax1)
            cbar.set_label('Temperature (K)')
            cbar.remove()
            ax1.legend()
            ax1.grid()
            ax1.set_xlim(ax.get_xlim()[0],ax.get_xlim()[1])
            ax1.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1])

        if save != None:
            plt.savefig(save,dpi=120,bbox_inches='tight')

        if plot_atmos:
            return fig, [ax,ax1]
        else:
            return fig, ax


    def generate_tex(self,filename = None, plot_atmos = False):

        if filename == None:
            filename = self.species

        citations = self.best_data_citations
        table = []
        for i in range(len(citations['citation'])):
            if citations['bibtex'][i]=='Resource not found.':
                bib_cite = citations['citation'][i]
            elif citations['bibtex'][i]=='No data':
                bib_cite = 'No data'
            else:
                ind1 = citations['bibtex'][i].find('{')+1
                ind2 = citations['bibtex'][i][ind1:].find(',')+ind1
                bib_cite = '\cite{'+citations['bibtex'][i][ind1:ind2]+'}'
            table.append([bib_cite,'%.1f'%citations['wavelength range'][i][0]+' - '+'%.1f'%citations['wavelength range'][i][1]])

        middle = tabulate(table, ['citation','wavelength range (nm)'], tablefmt="latex_raw")

        middle = r'\centering'+'\n'+middle
        middle = r'\begin{table}'+'\n'+middle
        middle = middle+'\n'+r'\end{table}'+'\n'

        fig,ax = self.plot(plot_atmos = plot_atmos,save = self.species+'.png')
        plt.close()

        figure = ''.join([r'\begin{figure}'+'\n',
                        r'\centering'+'\n',
                        r'\includegraphics[width=\textwidth]{'+self.species+'.png}\n',
                        r'\caption{'+self.species+'}\n',
                        r'\label{fig:'+self.species+'}\n',
                        r'\end{figure}'+'\n'])

        header = r'''\documentclass[11pt]{article}
\usepackage{geometry}
\geometry{margin=1in}
\usepackage{graphicx}
\nonstopmode
\begin{document}

'''
        footer = r'''

\bibliographystyle{unsrt}
'''+ r'\bibliography{'+self.species+'}'+'\n'+r'\end{document}'


        file = header+middle+'\n\n\n'+figure+footer

        fil = open(filename+'.tex','w')
        fil.write(header)
        fil.write(middle)
        fil.write('\n\n\n')
        fil.write(figure)
        fil.write(footer)
        fil.close()

        unique_bibs = list(set(citations['bibtex']))
        fil= open(filename+'.bib','w')
        for unique in unique_bibs:
            fil.write('\n')
            fil.write(unique)
            fil.write('\n')
        fil.close()


    def generate_pdf(self):
        cmds = ['pdflatex '+self.species+'.tex',
                'bibtex '+self.species+'.aux',
                'pdflatex '+self.species+'.tex',
                'pdflatex '+self.species+'.tex',
                'rm *.log *.blg *.bbl *.aux']
        for cmd in cmds:
            subprocess.call(cmd.split())
        subprocess.call('rm *.log *.blg *.bbl *.aux',shell=True)

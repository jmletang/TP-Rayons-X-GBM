import numpy as np
import subprocess
import requests
import os
import sys
new_path = r'/usr/local/lib'
sys.path.append(new_path)
new_path = r'/usr/local/lib/python3.8/site-packages'
sys.path.append(new_path)
new_path = r'/usr/local/progs/gate/gate-9.1-install/bin/'
sys.path.append(new_path)

import xraylib as xrl
import matplotlib.pyplot as plt
import matplotlib as mpl
import spekpy as sp
#from xpecgen import xpecgen as xg
#from tabletext import to_text

import gatetools.phsp as phsp
import gatetools as gt
import logging
#import pyvista
#
#def DisplayGate(material1,material2,cwd):
#    pl = pyvista.Plotter()
#    pl.import_vrml(f'{cwd}/g4_01.wrl')
#    pl.add_text(material1, font_size=30)
#    pl.show()
#
#    if os.path.exists(f'{cwd}/g4_03.wrl'):
#        pl = pyvista.Plotter()
#        pl.import_vrml(f'{cwd}/g4_03.wrl')
#        pl.add_text(material2, font_size=30)
#        pl.show()



def DisplayPhSp(material1,material2,cwd):
    logger = logging.getLogger(__name__)
 
    phsp_filename1 = f"{cwd}/phsp{material1}.root"
    phsp_filename2 = f"{cwd}/phsp{material2}.root"
    data1, read_keys1, m1 = phsp.load(phsp_filename1)
    data2, read_keys2, m2 = phsp.load(phsp_filename2)

    plt.close(1)
    fig = plt.figure(num=1,dpi=150,clear=True)
    mpl.rcParams.update({'font.size': 6})
    axs1 = plt.subplot(211)
    axs2 = plt.subplot(212)
    x1 = data1[:, read_keys1.index('Ekine')] * 1e3
    x2 = data2[:, read_keys1.index('Ekine')] * 1e3
    axs1.hist(x1, range=(0,x1.max()),bins=100,
                  histtype='stepfilled',
                  alpha=0.5)
    axs1.set_title(material1)
    axs1.set_xlabel('Énergie (keV)')
    axs1.set_xlim(xmin=0)
    axs1.set_ylabel('Coups')
    axs2.hist(x2, range=(0,x2.max()),bins=100,
                  histtype='stepfilled',
                  alpha=0.5)
    axs2.set_title(material2)
    axs2.set_xlabel('Énergie (keV)')
    axs2.set_xlim(xmin=0)
    axs2.set_ylabel('Coups')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    
def formula(compound):
    elt_CP = [xrl.AtomicNumberToSymbol(compound['Elements'][i]) for i in range(compound['nElements'])]
    atf_CP = [compound['massFractions'][i]/xrl.AtomicWeight(compound['Elements'][i])
              for i in range(compound['nElements'])]
    atfn_CP = [str(atf/min(atf_CP)) for atf in atf_CP]
    return "".join(np.core.defchararray.add(elt_CP,atfn_CP))

def GetDensity(material):
    if material=='H2C':
        cpH2C = xrl.GetCompoundDataNISTByName('Polyethylene')
        density = cpH2C['density']
    elif material=='H2O':
        cpH2O = xrl.GetCompoundDataNISTByName('Water, Liquid')
        density = cpH2O['density']
    else:
        Z=xrl.SymbolToAtomicNumber(material)
        density = xrl.ElementDensity(Z)
    return density

def simulator(N0=100, energy1=100, energy2=100, material1="H2C", thickness1=39, material2="Pb", thickness2=0.1):
    nmat=2
    N0limit = 10000
    if energy2==-1:
        energy2 = energy1
    #for i in range(2*nmat):
    #    if os.path.exists(f"attenuation/g4_0{i}.wrl"):
    #        os.remove(f"attenuation/g4_0{i}.wrl")
    Ndt=np.zeros((nmat))
    for i,energy,material,thickness in zip(range(nmat),[energy1,energy2],[material1,material2],[thickness1,thickness2]):
        if N0<=100:
            # Small simulation, we can run two Monte Carlo simulations with Gate to have some feedback
            with open("Gate.log", "w") as f:
              #subprocess.run(f'source ../.profile && Gate --qt -a [ENERGY,{energy}][THICKNESS,{thickness}][N0,{N0}][MATERIAL,{material}] mac/main.mac',
              subprocess.run(f'Gate --qt -a [ENERGY,{energy}][THICKNESS,{thickness}][N0,{N0}][MATERIAL,{material}] mac/main.mac',
                             shell=True, executable='/bin/bash',
                             cwd='attenuation',
                             stdout=f,
                             stderr=f,
                             universal_newlines=True,
                             bufsize=0)
            Ndt[i]=np.loadtxt('attenuation/output/fluence.txt')
        elif N0<=N0limit:
            # Large simulation, we can run two Monte Carlo simulations with disabled visu with Gate to have some feedback
            with open("Gate.log", "w") as f:
              #rm phsp 
              #subprocess.run(f'source ../.profile && Gate -a [ENERGY,{energy}][THICKNESS,{thickness}][N0,{N0}][MATERIAL,{material}] mac/main_novisu.mac',
              subprocess.run(f'cd attenuation && Gate -a [ENERGY,{energy}][THICKNESS,{thickness}][N0,{N0}][MATERIAL,{material}] mac/main_novisu.mac',
                             shell=True, executable='/bin/bash',
                             #cwd='attenuation',
                             stdout=f,
                             stderr=f,
                             universal_newlines=True,
                             bufsize=0)
            Ndt[i]=np.loadtxt('attenuation/output/fluence.txt')
        else:
            density = GetDensity(material)
            Ndt[i] = np.random.poisson(N0*np.exp(-thickness*density*xrl.CS_Total_CP(material, energy)))            
        print(f'{int(Ndt[i])} photons ont traversé la plaque de {thickness} cm de {material} sur {N0} envoyés')
    if N0>100 and N0<=N0limit:
        DisplayPhSp(material1,material2,'attenuation/output')
    elif N0>N0limit:
        print(f'La visualisation est désactivée au delà de {N0limit} photons émis.')

def mu(material='H2C'):
    energy_range = np.arange(5.,800., 0.1, dtype=np.double)
    density = GetDensity(material)
    #print(f'density {material} = {density}')
    mu_rho = [xrl.CS_Total_CP(material, E) * density for E in energy_range]
    mu_rho_Photo = [xrl.CS_Photo_CP(material, E) * density for E in energy_range]
    mu_rho_Compt = [xrl.CS_Compt_CP(material, E) * density for E in energy_range]
    mu_rho_Rayl = [xrl.CS_Rayl_CP(material, E) * density for E in energy_range]
    plt.close(1)
    fig = plt.figure(num=1,dpi=150,clear=True)
    mpl.rcParams.update({'font.size': 6})
    axMW = plt.subplot(111)
    axMW.plot(energy_range, mu_rho,color="black",linewidth=2.,linestyle="-",label='Total')
    axMW.plot(energy_range, mu_rho_Photo,color="red",linewidth=2.,linestyle="-",label='Photoelectric')
    axMW.plot(energy_range, mu_rho_Compt,color="blue",linewidth=2.,linestyle="-",label='Compton')
    axMW.plot(energy_range, mu_rho_Rayl,color="green",linewidth=2.,linestyle="-",label='Rayleigh')
    axMW.set_xscale('log')
    axMW.set_yscale('log')
    axMW.set_xlim(np.min(energy_range),np.max(energy_range))
    axMW.set_ylim(1e-2,1e4)
    plt.legend(loc='center right', frameon=True)
    plt.xlabel('Énergie (keV)')
    plt.ylabel("Coefficient d'atténuation linéique (cm$^{-1}$)")
    axMW.grid(which='major', axis='x', linewidth=0.5, linestyle='-', color='0.75')
    axMW.grid(which='minor', axis='x', linewidth=0.3, linestyle='-', color='0.75')
    axMW.grid(which='major', axis='y', linewidth=0.5, linestyle='-', color='0.75')
    axMW.grid(which='minor', axis='y', linewidth=0.3, linestyle='-', color='0.75')
    axMW.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%d"))
    #axMW.xaxis.set_minor_formatter(mpl.ticker.FormatStrFormatter("%d"))
    axMW.grid(True)
    #symbol=xrl.AtomicNumberToSymbol(material)
    axMW.set_title("%s" % material, va='bottom')
    #plt.savefig('mu_over_rho_W.pdf', format='PDF')
    text=axMW.text(np.min(energy_range),1e4, "", va="top", ha="left")
    def onclick(event):
        energy = np.round(event.xdata*10)*0.1
        energyidx = int(np.where(np.min(np.abs(energy_range-energy))==np.abs(energy_range-energy))[0])
        tx = 'Le coefficient d\'atténuation linéique de ' + material + ' à %.1f keV est %1.4e cm$^{-1}$\n(Rayleigh %1.4e cm$^{-1}$, Photoelectric %1.4e cm$^{-1}$, Compton %1.4e cm$^{-1}$)'%(energy,mu_rho[energyidx],mu_rho_Rayl[energyidx],mu_rho_Photo[energyidx],mu_rho_Compt[energyidx])
        text.set_text(tx)
        text.set_x(axMW.get_xlim()[0])
        text.set_y(axMW.get_ylim()[1])
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

def spectrum(E0,Mat_Z,Mat_X):
    xrs=sp.Spek(kvp=E0,th=12,dk=0.5,physics="casim",mu_data_source="pene",mas=1.0)
    #xrs=xg.calculate_spectrum(E0,12,3,100,epsrel=0.5,monitor=None,z=74)
    #Inherent filtration: 1.2mm Al + 100cm Air
    xrs.filter('Al',1.2).filter('Air',1000)
    #mu_Al=xg.get_mu(13)
    #mu_Cu=xg.get_mu(29)
    #mu_Pb=xg.get_mu(82)
    #xrs.attenuate(0.12,mu_Al)
    #xrs.attenuate(100,xg.get_mu("air"))

    #fluence_to_dose=xg.get_fluence_to_dose()
    #xrs.set_norm(value=1.0,weight=fluence_to_dose)
    #Attenuation
    #if Mat_Z>0: #Atomic number
    #    dMat = xrl.ElementDensity(Mat_Z)
    #    fMat = xrl.AtomicNumberToSymbol(Mat_Z)
    #    xrs.attenuate(0.1*Mat_X,xg.get_mu(Mat_Z))
    #else: #-1 == 'Water, Liquid'    
    #    mH2O = 2. * xrl.AtomicWeight(1) + xrl.AtomicWeight(8)
    #    wH = 0.1 * Mat_X * 2. * xrl.AtomicWeight(1) / (xrl.ElementDensity(1) * mH2O)
    #    wO = 0.1 * Mat_X * xrl.AtomicWeight(8) / (xrl.ElementDensity(8) * mH2O)
    #    xrs.attenuate(wH,xg.get_mu(1))
    #    xrs.attenuate(wO,xg.get_mu(8))
    xrs.filter(Mat_Z,Mat_X)

    #Get the figures
    Nr_Photons = "%.4g photons" % (np.sum(xrs.get_spk()))
    Average_Energy = "%.2f keV" % (xrs.get_emean())
    Dose = "%.3g mGy (air kerma)" % (xrs.get_kerma())
    HVL_Al_text = "%.3f mm(%s)" % (xrs.get_hvl1(matl=Mat_Z),Mat_Z)
    print(f"Dose = {Dose} per [mA·s] @ 1m")
    print(f"Fluence = {Nr_Photons} per [cm²·mA·s] @ 1m")
    print(f"Average energy = {Average_Energy}")
    print(f"Half-value Layer (for dose) = {HVL_Al_text}")
    #Nr_Photons = "%.3g" % (xrs.get_norm())
    #Average_Energy = "%.1f keV" % (xrs.get_norm(lambda x:x)/xrs.get_norm())
    #Dose = "%.2g mGy" % (xrs.get_norm(fluence_to_dose))
    #HVL_Al=xrs.hvl(0.5,fluence_to_dose,mu_Al)
    #HVL_Cu=xrs.hvl(0.5,fluence_to_dose,mu_Cu)
    #HVL_Pb=xrs.hvl(0.5,fluence_to_dose,mu_Pb)
    #HVL_text = "%.3f mm (Al) ou %.3f mm (Cu) ou %.3f mm (Pb)" % (10*HVL_Al,10*HVL_Cu,10*HVL_Pb)
    #a = [["Dose à 1m", Dose],["Nombre total de photons", Nr_Photons],
    #     ["Énergie moyenne des photons",Average_Energy],
    #     ["Couche de Demi-Atténuation", HVL_text]]
    #print(to_text(a))
    #print(a)
    
    #(x2,y2) = xrs.get_points()
    (x2,y2) = xrs.get_spectrum()

    plt.close(2)
    plt.figure(num=2,dpi=150,clear=True)
    mpl.rcParams.update({'font.size': 6})
    axMW = plt.subplot(111)
    axMW.plot(x2,y2)
    axMW.set_xlim(3,E0)
    axMW.set_ylim(0,)
    plt.xlabel("Énergie [keV]")
    plt.ylabel("Nombre de photons par [keV·cm²·mGy] @ 1m")
    axMW.grid(which='major', axis='x', linewidth=0.5, linestyle='-', color='0.75')
    axMW.grid(which='minor', axis='x', linewidth=0.2, linestyle='-', color='0.85')
    axMW.grid(which='major', axis='y', linewidth=0.5, linestyle='-', color='0.75')
    axMW.grid(which='minor', axis='y', linewidth=0.2, linestyle='-', color='0.85')
    axMW.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%d"))
    axMW.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2g"))
    axMW.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
    axMW.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
    axMW.grid(True)
    plt.show()

def simulator_scatter(N0=100, position=39, sdd=1000):
    if os.path.exists(f"diffuse/g4_00.wrl"):
        os.remove(f"diffuse/g4_00.wrl")
    if os.path.exists(f"diffuse/g4_01.wrl"):
        os.remove(f"diffuse/g4_01.wrl")
    if os.path.exists(f"diffuse/output/fluence.mhd"):
        os.remove(f"diffuse/output/fluence.mhd")
    if os.path.exists(f"diffuse/output/fluence.raw"):
        os.remove(f"diffuse/output/fluence.raw")
    if N0<=100:
        macro='mac/main_visu.mac'
    else:
        macro='mac/main_novisu.mac'
    # Small simulation, we can run two Monte Carlo simulations with Gate to have some feedback
    with open("Gate.log", "w") as f:
      #subprocess.run(f'source ../.profile && Gate {"--qt" if N0<=100 else ""} -a [N0,{N0}][PLACEMENT,{position}][SDD,{sdd}] {macro}',
      subprocess.run(f'Gate {"--qt" if N0<=100 else ""} -a [N0,{N0}][PLACEMENT,{position}][SDD,{sdd}] {macro}',
                     shell=True, executable='/bin/bash',
                     cwd='diffuse',
                     stdout=f,
                     stderr=f,
                     universal_newlines=True,
                     bufsize=0)
    if N0>100:
        print('La visualisation 3D est désactivée au delà de 100 photons émis.')

    fluence = np.fromfile('diffuse/output/fluence.raw',dtype=np.float32)
    fluence = fluence.reshape((1536,2048))
    binning = 16
    shape = (1536//binning, binning, 2048//binning, binning)
    fluence = fluence.reshape(shape).sum(-1).sum(1)
    
    plt.close(12)
    plt.figure(num=12,figsize=(6,2), dpi=150)
    plt.subplot(121)
    plt.imshow(fluence, cmap=plt.cm.gray,extent=(0,2047,0,1535))
    plt.title(f'À {position} cm du détecteur (binning {binning}x{binning})', fontsize=8)
    plt.subplot(122)
    fluence = fluence[512//binning:1024//binning,:]
    plt.plot(np.arange(0,2048,binning),np.mean(fluence,axis=0))
    plt.title(f'Profil horizontal moyen', fontsize=8)
    plt.show()

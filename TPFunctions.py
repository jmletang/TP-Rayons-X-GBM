import numpy as np
import subprocess
import requests
import os
from IPython.display import IFrame
from IPython.core.display import display, HTML
from bs4 import BeautifulSoup
import xraylib as xrl
import matplotlib.pyplot as plt
import matplotlib as mpl
from xpecgen import xpecgen as xg
from tabletext import to_text

def DisplayGate(material1,material2):
    pxheight=800
    # Create html file for the simulation
    html_file = open("simulation.html", "w")
    with open(f'attenuation/g4_01.wrl', 'r') as file:
        wrl_data = file.read()
    data = {'input_encoding': 'CLASSIC',
            'output_encoding': 'HTML5',
            'input_code': wrl_data}
    API_ENDPOINT = "https://doc.instantreality.org/tools/x3d_encoding_converter/convert/"
    r = requests.post(url=API_ENDPOINT, data=data)
    soup = BeautifulSoup(r.text, features='html.parser')
    soup = BeautifulSoup(soup.find('div', {'class': 'source'}).text, features='html.parser')
    soup.link["href"]="https://www.x3dom.org/x3dom/release/x3dom.css"
    soup.script["src"]="https://www.x3dom.org/x3dom/release/x3dom.js"
    soup.x3d["width"]="100%"
    soup.x3d["height"]=f"{pxheight}px"
    soup.style="float:left;"
    if os.path.exists('attenuation/g4_03.wrl'):
        with open('attenuation/g4_03.wrl', 'r') as file:
            wrl_data = file.read()
        data = {'input_encoding': 'CLASSIC',
                'output_encoding': 'HTML5',
                'input_code': wrl_data}
        API_ENDPOINT = "https://doc.instantreality.org/tools/x3d_encoding_converter/convert/"
        r = requests.post(url=API_ENDPOINT, data=data)
        soup2 = BeautifulSoup(r.text, features='html.parser')
        soup2 = BeautifulSoup(soup2.find('div', {'class': 'source'}).text, features='html.parser')
        soup2.x3d["width"]="100%"
        soup2.x3d["height"]=f"{pxheight}px"

    style = soup.new_tag("style")
    style.string = ".split { height: 100%; width: 50%; position: fixed; top: 0; } .left { left: 0; } .right { right: 0; border-left: 3px solid black; }"
    soup.head.append(style)

    left = soup.new_tag("div")
    left["class"] = "split left"
    left.append(material1)
    left.append(soup.x3d)
    right = soup.new_tag("div")
    right["class"] = "split right"
    right.append(material2)
    right.append(soup2.x3d)
    soup.body.clear()
    soup.body.append(left)
    soup.body.append(right)

    html_file.write(soup.prettify())
    html_file.close()
    display(HTML("<style>.container { width:100% !important; }</style>"))
    display(IFrame(src='./simulation.html', width="100%", align="left", height=pxheight))

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
    else:
        Z=xrl.SymbolToAtomicNumber(material)
        density = xrl.ElementDensity(Z)
    return density

def simulator(N0=100, energy1=100, energy2=100, material1="H2C", thickness1=39, material2="Pb", thickness2=0.1):
    nmat=2
    for i in range(2*nmat):
        if os.path.exists(f"attenuation/g4_0{i}.wrl"):
            os.remove(f"attenuation/g4_0{i}.wrl")
    Ndt=np.zeros((nmat))
    for i,energy,material,thickness in zip(range(nmat),[energy1,energy2],[material1,material2],[thickness1,thickness2]):
        if N0<=100:
            # Small simulation, we can run two Monte Carlo simulations with Gate to have some feedback
            with open("Gate.log", "w") as f:
              subprocess.run(f'source ~/.profile && Gate -a [ENERGY,{energy}][THICKNESS,{thickness}][N0,{N0}][MATERIAL,{material}] mac/main.mac',
                             shell=True, executable='/bin/bash',
                             cwd='attenuation',
                             stdout=f,
                             stderr=f,
                             universal_newlines=True,
                             bufsize=0)
            Ndt[i]=np.loadtxt('attenuation/output/fluence.txt')
        else:
            density = GetDensity(material)
            Ndt[i] = np.random.poisson(N0*np.exp(-thickness*0.1*density*xrl.CS_Total_CP(material, energy)))
        print(f'{int(Ndt[i])} photons ont traversé la plaque de {thickness} cm de {material} sur {N0} envoyés')
    if N0<=100:
        DisplayGate(material1,material2)
    else:
        print('La visualisation est désactivée au delà de 100 photons émis.')

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
    xrs=xg.calculate_spectrum(E0,12,3,100,epsrel=0.5,monitor=None,z=74)
    #Inherent filtration: 1.2mm Al + 100cm Air
    mu_Al=xg.get_mu(13)
    xrs.attenuate(0.12,mu_Al)
    xrs.attenuate(100,xg.get_mu("air"))
    fluence_to_dose=xg.get_fluence_to_dose()
    xrs.set_norm(value=0.146,weight=fluence_to_dose)
    #Attenuation
    if Mat_Z>0: #Atomic number
        dMat = xrl.ElementDensity(Mat_Z)
        fMat = xrl.AtomicNumberToSymbol(Mat_Z)
        xrs.attenuate(0.1*Mat_X,xg.get_mu(Mat_Z))
    else: #-1 == 'Water'
        mH2O = 2. * xrl.AtomicWeight(1) + xrl.AtomicWeight(8)
        wH = 0.1 * Mat_X * 2. * xrl.AtomicWeight(1) / (xrl.ElementDensity(1) * mH2O)
        wO = 0.1 * Mat_X * xrl.AtomicWeight(8) / (xrl.ElementDensity(8) * mH2O)
        xrs.attenuate(wH,xg.get_mu(1))
        xrs.attenuate(wO,xg.get_mu(8))
    #Get the figures
    Nr_Photons = "%.4g" % (xrs.get_norm())
    Average_Energy = "%.2f keV" % (xrs.get_norm(lambda x:x)/xrs.get_norm())
    Dose = "%.3g mGy" % (xrs.get_norm(fluence_to_dose))
    HVL_Al=xrs.hvl(0.5,fluence_to_dose,mu_Al)
    HVL_Al_text = "%.2f mm (Al)" % (10*HVL_Al)
    a = [["Dose à 1m", Dose],["Nombre total de photons", Nr_Photons],
         ["Énergie moyenne des photons",Average_Energy],["Couche de Demi-Atténuation", HVL_Al_text]]
    print(to_text(a))
    (x2,y2) = xrs.get_points()
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
    axMW.xaxis.set_minor_locator(AutoMinorLocator())
    axMW.yaxis.set_minor_locator(AutoMinorLocator())
    axMW.grid(True)
    plt.show()


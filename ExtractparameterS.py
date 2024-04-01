import re
from io import StringIO
import numpy as np
import  matplotlib.pyplot as plt

file_name="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wf1_Nf16_vbg_1V_Spar_HF_cold_0V.mdm"
id_vd="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wfp5_Nf32_vbg_0V_idvd_0p8V.mdm"
extract_file="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wfp5_Nf16_BGF_vbg_0V_Spar_HF_cold_0V.mdm"#vd=0,vbg=0, vg =-0.2 ->0.9
id_vg="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wfp5_Nf32_vbg_0V_idvg_50m.mdm"
def mdm_reader_index(fname,measurement_index):
    f = open(fname, "r")    
    text = f.read()
    ldata = []
    i=0
    for begin_match, end_match in zip(re.finditer("BEGIN_DB", text), re.finditer("END_DB", text)):
        
        s = begin_match.end()
        e = end_match.start()
        body = text[s:e]
        
        if i==measurement_index:
            data,info = read_block(body)
            ldata.append(data)
        i+=1
        # convert to numpy array
    adata = np.array(ldata)
    f.close()
    return adata,info
def mdm_reader_extract(fname):
    f = open(fname, "r")    
    text = f.read()
    ldata = []
    i=0
    for begin_match, end_match in zip(re.finditer("BEGIN_DB", text), re.finditer("END_DB", text)):
        
        s = begin_match.end()
        e = end_match.start()
        body = text[s:e]
        
       
        data,info = read_block(body)
        ldata.append(data)
        
        # convert to numpy array
    adata = np.array(ldata)
    f.close()
    return adata
def read_block(body):
    lines = body.strip().split("\n")
    output = StringIO()
    info=[]
    for line in lines:
        tline= line.strip()
        if tline.startswith("ICCAP_VAR") or tline.startswith("#"):
            info.append(tline)
            new_line = "%"+ tline
        else:
            
            new_line = tline
        # only write non blank lines
        if (new_line):
            output.write(new_line + '\n')
    output.seek(0)
    data = np.loadtxt(output, comments='%')
    return data,info




#shape of read: triple array: premiere dim= les differentes mesure (par exemple pour different vg), 
# la deuxieme= array contenant les mesures sur le parametre que l'on sweep
#la troisieme= les differentes mesures pour une valeur de sweep fixée, par ex id,ig, S11, S12 ect


def plot_id_vd(file_name, measurement_index, id_index , vd_index):
    data,info=mdm_reader_index(file_name,measurement_index)
    measurement=data[0]
    
   
    vg= info[0].split()[-1]
    
    vd = [column[ vd_index] for column in measurement]
    id = [column[id_index] for column in measurement]

    # Tracer les valeurs de la deuxième colonne par rapport à celles de la première colonne
    plt.plot(vd, id , linestyle='-')
    plt.xlabel('vd[V]')
    plt.ylabel('id[V]')
    plt.title('id vs vd curve for vg={} [V]'.format(vg))
    plt.grid(True)
    plt.show()
#plot_id_vd(id_vd,0,1,0)
def plot_id_vg( file_name,measurement_index, id_index , vg_index):
    data,info=mdm_reader_index(file_name,measurement_index)
    measurement=data[0]
    
   
    vd= info[0].split()[-1]
    
    vg = [column[vg_index] for column in measurement]
    id = [column[id_index] for column in measurement]

    # Tracer les valeurs de la deuxième colonne par rapport à celles de la première colonne
    plt.plot(vg, id , linestyle='-')
    plt.xlabel('vg[V]')
    plt.ylabel('id[V]')
    plt.title('id vs vd curve for vd={} [V]'.format(vd))
    plt.grid(True)
    plt.show()
    deriv_1=np.gradient(id,0.01)
    deriv_2=np.gradient(deriv_1,0.01)
    indice_max=np.argmax(deriv_2)
    print("Vth=",vg[indice_max])

#plot_id_vg(id_vg,0,1,0)
def plot_S_parameter( measurement_index):
    data,info=_index(file_name,measurement_index)
    measurement=data[0]
    
    vd= info[0].split()[-1]
    vg= info[1].split()[-1]
    

    freq = [column[0] for column in measurement]
    S11 = 20*np.log10(np.sqrt(np.square(measurement[:,3])+np.square(measurement[:,4])))
    S12 = 20*np.log10(np.sqrt(np.square(measurement[:,5])+np.square(measurement[:,6])))
    S21 = 20*np.log10(np.sqrt(np.square(measurement[:,7])+np.square(measurement[:,8])))
    S22 = 20*np.log10(np.sqrt(np.square(measurement[:,9])+np.square(measurement[:,10])))
    
    plt.plot(freq, S11 , linestyle='-',label="S11, Vd={},Vg={}".format(vd,vg))
    plt.plot(freq, S21 , linestyle='-',label="S21,,Vd={},Vg={}".format(vd,vg))
    plt.xlabel('freq[Hz]')
    plt.ylabel('magnitude[dB]')
    plt.title('S parameter')
    plt.text(1, 20, info[0], fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
    plt.legend()
    plt.show()
    
    
#plot_S_parameter(10)
#Sij ici est nombre imaginaire, pas amplitude
def calcul_Y11(S11, S12, S21, S22):
    return ((1 + S22) * (1 - S11) + S12 * S21) / ((1 + S22) * (1 + S11) - S12 * S21)
def calcul_Y12(S11, S12, S21, S22):
    return -2*S12/((1 + S22) * (1 + S11) - S12 * S21)
def calcul_Y21(S11, S12, S21, S22):
    return -2*S21/((1 + S22) * (1 + S11) - S12 * S21)
def calcul_Y22(S11, S12, S21, S22):
    return ((1 + S11) * (1 - S22) + S12 * S21) / ((1 + S22) * (1 + S11) - S12 * S21)



def plot_Y_parameter( measurement_index):
    data,info=mdm_reader_index(file_name,measurement_index)
    measurement=data[0]
    
    vd= info[0].split()[-1]
    vg= info[1].split()[-1]
    

    freq = [column[0] for column in measurement]
    S11 = np.array(measurement[:,3]+1j*measurement[:,4],dtype="complex")
    S12 = np.array(measurement[:,5]+1j*measurement[:,6],dtype="complex")
    S21 = np.array(measurement[:,7]+1j*measurement[:,8],dtype="complex")
    S22 = np.array(measurement[:,9]+1j*measurement[:,10],dtype="complex")
    
    
    Y11= calcul_Y11(S11,S12,S21,S22)
    Y21= calcul_Y21(S11,S12,S21,S22)
    Y12= calcul_Y12(S11,S12,S21,S22)
    Y22= calcul_Y22(S11,S12,S21,S22)
    
    
    plt.plot(freq, Y11.real , linestyle='-',label="Y11, Vd={},Vg={}".format(vd,vg))
    plt.plot(freq, Y21.real , linestyle='-',label="Y21,,Vd={},Vg={}".format(vd,vg))
    plt.xlabel('freq[Hz]')
    plt.ylabel('amplitude')
    plt.title('Y parameter')
    plt.text(1, 20, info[0], fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
    plt.legend()
    plt.show()
    
    
#plot_Y_parameter(10)


def calcul_Z11(S11, S12, S21, S22):
    return ((1 + S11) * (1 - S22) + S12 * S21) / ((1 -S11) * (1 - S22) - S12 * S21)
def calcul_Z12(S11, S12, S21, S22):
    return 2*S12/((1 - S11) * (1 -S22) - S12 * S21)
def calcul_Z21(S11, S12, S21, S22):
   return 2*S21/((1 - S11) * (1 -S22) - S12 * S21)
def calcul_Z22(S11, S12, S21, S22):
    return  ((1 + S22) * (1 - S11) + S12 * S21) / ((1 -S11) * (1 - S22) - S12 * S21)


def plot_Z_parameter( measurement_index):
    data,info=mdm_reader_index(file_name,measurement_index)
    measurement=data[0]
    
    vd= info[0].split()[-1]
    vg= info[1].split()[-1]
    

    freq = [column[0] for column in measurement]
    S11 = np.array(measurement[:,3]+1j*measurement[:,4],dtype="complex")
    S12 = np.array(measurement[:,5]+1j*measurement[:,6],dtype="complex")
    S21 = np.array(measurement[:,7]+1j*measurement[:,8],dtype="complex")
    S22 = np.array(measurement[:,9]+1j*measurement[:,10],dtype="complex")
    
    
    Z11= calcul_Z11(S11,S12,S21,S22)
    Z21= calcul_Z21(S11,S12,S21,S22)
    Z12= calcul_Z12(S11,S12,S21,S22)
    Z22= calcul_Z22(S11,S12,S21,S22)
    
    
    plt.plot(freq, Z11.real , linestyle='-',label="Z11, Vd={},Vg={}".format(vd,vg))
    plt.plot(freq, Z21.real , linestyle='-',label="Z21,,Vd={},Vg={}".format(vd,vg))
    plt.xlabel('freq[Hz]')
    plt.ylabel('amplitude')
    plt.title('Z parameter')
    plt.text(1, 20, info[0], fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
    plt.legend()
    plt.show()
    
    
#plot_Z_parameter(10)

def extract_Zparameter_Vg(file_name):
    data=mdm_reader_extract(file_name)
    
    vth=0.29
    measurement=data[:,0]
    print(measurement)
    measurement=measurement[4:]#on Retire les vg<vth
    
    S11 = np.array(measurement[:,3]+1j*measurement[:,4],dtype="complex")
    S12 = np.array(measurement[:,5]+1j*measurement[:,6],dtype="complex")
    S21 = np.array(measurement[:,7]+1j*measurement[:,8],dtype="complex")
    S22 = np.array(measurement[:,9]+1j*measurement[:,10],dtype="complex")
    

    Z11= calcul_Z11(S11,S12,S21,S22)
    Z21= calcul_Z21(S11,S12,S21,S22)
    Z12= calcul_Z12(S11,S12,S21,S22)
    Z22= calcul_Z22(S11,S12,S21,S22)
    print(Z22.shape)
    
    vg=[ 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
    vec_th=np.ones(len(vg))*vth
    plt.plot(1/(vg-vec_th), Z11.real , linestyle='-',label="Z22real")
    plt.grid(True)
    plt.xlabel('1/Vg-vth')
    plt.ylabel('Real Z22')
    plt.title('Z parameter extraction of gdsin ')
    plt.legend()
    plt.show()
    
    
  
extract_Zparameter_Vg(extract_file)
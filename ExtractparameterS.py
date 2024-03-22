import re
from io import StringIO
import numpy as np
import  matplotlib.pyplot as plt

file_name="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wf1_Nf16_vbg_1V_Spar_HF_cold_0V.mdm"
def mdm_reader(fname):
    f = open(fname, "r")    
    text = f.read()
    ldata = []
    for begin_match, end_match in zip(re.finditer("BEGIN_DB", text), re.finditer("END_DB", text)):
        s = begin_match.end()
        e = end_match.start()
        body = text[s:e]
        data = read_block(body)
        ldata.append(data)
        # convert to numpy array
    adata = np.array(ldata)
    f.close()
    return adata
def read_block(body):
    lines = body.strip().split("\n")
    output = StringIO()
    for line in lines:
        tline= line.strip()
        if tline.startswith("ICCAP_VAR") or tline.startswith("#"):
            new_line = "%"+ tline
        else:
            new_line = tline
        # only write non blank lines
        if (new_line):
            output.write(new_line + '\n')
    output.seek(0)
    data = np.loadtxt(output, comments='%')
    return data


read=mdm_reader(file_name)
#shape of read: triple array: premiere dim= les differentes mesure (par exemple pour different vg), 
# la deuxieme= array contenant les mesures sur le parametre que l'on sweep
#la troisieme= les differentes mesures pour une valeur de sweep fixée, par ex id,ig, S11, S12 ect


def plot_i_v(readed_file, measurement_index, id_index , vd_index):
    measurement=read[measurement_index]
    
    vd = [column[ vd_index] for column in measurement]
    id = [column[id_index] for column in measurement]

    # Tracer les valeurs de la deuxième colonne par rapport à celles de la première colonne
    plt.plot(vd, id , linestyle='-')
    plt.xlabel('vd[V]')
    plt.ylabel('id[V]')
    plt.title('id vs vd curve for vg=0.8 [V]')
    plt.grid(True)
    plt.show()


def plot_S_parameter(readed_file, measurement_index):
    measurement=read[measurement_index]
    freq = [column[0] for column in measurement]
    S11 = 20*np.log10(np.sqrt(np.square(measurement[:,3])+np.square(measurement[:,4])))
    S12 = 20*np.log10(np.sqrt(np.square(measurement[:,5])+np.square(measurement[:,6])))
    S21 = 20*np.log10(np.sqrt(np.square(measurement[:,7])+np.square(measurement[:,8])))
    S22 = 20*np.log10(np.sqrt(np.square(measurement[:,9])+np.square(measurement[:,10])))
    
    plt.plot(freq, S11 , linestyle='-',label="S11")
    plt.plot(freq, S21 , linestyle='-',label="S21")
    plt.xlabel('freq[Hz]')
    plt.ylabel('magnitude[dB]')
    plt.title('S parameter')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    
plot_S_parameter(read, 1)
    
    
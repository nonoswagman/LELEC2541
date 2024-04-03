import re
from io import StringIO
import numpy as np
import  matplotlib.pyplot as plt
from sympy import symbols, solve

id_vd="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wfp5_Nf32_vbg_0V_idvd_0p8V.mdm"
extract_file="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wfp5_Nf32_vbg_0V_Spar_HF_cold_0V.mdm"#vd=0,vbg=0, vg =-0.2 ->0.9
open_measure="1M/data_MPW2230_Die1_20200708/Die1_SLVTL20_Wfp5_Nf32_OPEN_Spar_HF.mdm"
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
    plt.title('id vs vg curve for vd={} [V]'.format(vd))
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
    # return [[Y11,Y12],[Y21,Y22]]
    
    
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
def correct_Z_parameter_for_extraction(Y11,Y12,Y21,Y22,freq_indice,open_file,len_vg):
    data=mdm_reader_extract(open_file)
    measurement=data[:,freq_indice]
    #Open mesure fait pour un seul vg=0 ->on retire le meme open a toute les valeurs de Y11 pour les differents Vg
    S11 = np.array(measurement[:,3]+1j*measurement[:,4],dtype="complex")
    S12 = np.array(measurement[:,5]+1j*measurement[:,6],dtype="complex")
    S21 = np.array(measurement[:,7]+1j*measurement[:,8],dtype="complex")
    S22 = np.array(measurement[:,9]+1j*measurement[:,10],dtype="complex")
    

    Y11_op= calcul_Y11(S11,S12,S21,S22)
    Y21_op= calcul_Y21(S11,S12,S21,S22)
    Y12_op= calcul_Y12(S11,S12,S21,S22)
    Y22_op= calcul_Y22(S11,S12,S21,S22)
   
    Y11_tot=Y11-Y11_op
    
    Y22_tot=Y22-Y22_op
    Y21_tot=Y21-Y21_op
    Y12_tot=Y12-Y12_op
    Y=np.array([[Y11_tot,Y12_tot],[Y21_tot,Y22_tot]])
    Z_matrix_tot=np.zeros((2, 2, len(Y11)), dtype=complex)
    for i in range(len(Y11)):
        sub_matrix=Y[:,:,i] # On prend la matrice Y pour un vg i fixé et on l'inverse puis la remet dans Z à la bonne place
        inv = np.linalg.inv(sub_matrix)
        Z_matrix_tot[:,:,i]=inv
    return Z_matrix_tot[0][0],Z_matrix_tot[0][1],Z_matrix_tot[1][0],Z_matrix_tot[1][1]
 
def extract_Extrinsic_resistance(file_name,freq_indice_range,open_file):
    'Retourne Rge,Lge,Rse,Lse,Rde,Lde'
    data=mdm_reader_extract(file_name)
    measurement1=data[:,freq_indice_range[0]]
    measurement2=data[:,freq_indice_range[-1]]
    
    omega1=measurement1[0][0]
    omega2=measurement2[0][0]
    print("Les fréquence utilisées sur lesquelles on fait la moyenne des Z parametres sont :", (omega1,omega2) )
    
    vth=0.29
    vg=[ 0.75, 0.8, 0.85, 0.9]
    freq_list= list(range(freq_indice_range[0], freq_indice_range[1] + 1))
    Z11_mean=np.zeros(len(vg),dtype="complex")
    Z12_mean=np.zeros(len(vg),dtype="complex")
    Z21_mean=np.zeros(len(vg),dtype="complex")
    Z22_mean=np.zeros(len(vg),dtype="complex")
    for freq_indice in freq_list:
        
    
        measurement=data[:,freq_indice]
        measurement=measurement[(len(measurement)-len(vg)):]#on Retire les vg trop faible
        
        S11 = np.array(measurement[:,3]+1j*measurement[:,4],dtype="complex")
        S12 = np.array(measurement[:,5]+1j*measurement[:,6],dtype="complex")
        S21 = np.array(measurement[:,7]+1j*measurement[:,8],dtype="complex")
        S22 = np.array(measurement[:,9]+1j*measurement[:,10],dtype="complex")
        

        Y11= calcul_Y11(S11,S12,S21,S22)
        Y21= calcul_Y21(S11,S12,S21,S22)
        Y12= calcul_Y12(S11,S12,S21,S22)
        Y22= calcul_Y22(S11,S12,S21,S22)
        #! Yij et Zij sont des vecteur contenant les differentes valeurs de Yij/Zij poue different vf, a une freq fixée
        # On fait une moyenne pour retirer la dépendance en freq de Re(Ze)
        Z11,Z12,Z21,Z22=correct_Z_parameter_for_extraction(Y11,Y12,Y21,Y22,freq_indice,open_file,len(vg))
        Z11_mean+=Z11/len(freq_list)
        Z12_mean+=Z12/len(freq_list)
        Z21_mean+=Z21/len(freq_list)
        Z22_mean+=Z22/len(freq_list)
        
    ################################################Extraction partie réelle
    
    slope_11 = (Z11_mean[-1].real - Z11_mean[0].real) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
    
    extrapolated_Z11_real = Z11_mean[0].real + slope_11* (0 - 1/(vg[0]-vth))
    
    slope_12 = (Z12_mean[-1].real - Z12_mean[0].real) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
   
    extrapolated_Z12_real = Z12_mean[0].real + slope_12 * (0 - 1/(vg[0]-vth))
    
    slope_21 = (Z21_mean[-1].real - Z21_mean[0].real) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
    
    extrapolated_Z21_real = Z21_mean[0].real + slope_21 * (0 - 1/(vg[0]-vth))
    
    slope_22 = (Z22_mean[-1].real - Z22_mean[0].real) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
    
    extrapolated_Z22_real = Z22_mean[0].real + slope_22 * (0 - 1/(vg[0]-vth))
    
    
    
    
    
    
    Rg, Rs, Rd = symbols('Rg Rs Rd')

    # Déclaration des équations
    eq1 = Rs + Rg - extrapolated_Z11_real
    eq2 = Rs - extrapolated_Z12_real
    eq3 = Rd + Rs - extrapolated_Z22_real

    # Résolution du système d'équations
    solution_R = solve((eq1, eq2, eq3), (Rg, Rs, Rd))
    
    
    
    
    
    return solution_R[Rg],solution_R[Rs],solution_R[Rd]

    # print("La valeur extrapolée de Z11 pour x=0 est :", extrapolated_Z11)
    # vec_th=np.ones(len(vg))*vth
    # plt.plot(1/(vg-vec_th), Z11_mean.real , linestyle='-',label="Z22real")
    # plt.plot(0, extrapolated_Z11, marker='o', linestyle='--', color='red', label='Extrapolation')
    # plt.grid(True)
    # plt.xlabel('1/Vg-vth')
    # plt.ylabel('Real Z22 ')
    # plt.title('Z parameter extraction of gdsin after correction ')
    # plt.legend()
    # plt.show()
    
def extract_Extrinsic_Inductance(file_name,freq_indice,open_file):
    "ici comme on doit connaitre L indépendemment de f, on doit calculer pour une frequence fixée, pas une moyenne de f comme pour R"
    data=mdm_reader_extract(file_name)
    
    vth=0.29
    vg=[ 0.75, 0.8, 0.85, 0.9]
    measurement=data[:,freq_indice]
    
    omega=measurement[0][0]
    print("La fréquence utilisée pour le calcul des inductances séries extrinsèques est:",omega )
    measurement=measurement[(len(measurement)-len(vg)):]#on Retire les vg trop faible
    
    S11 = np.array(measurement[:,3]+1j*measurement[:,4],dtype="complex")
    S12 = np.array(measurement[:,5]+1j*measurement[:,6],dtype="complex")
    S21 = np.array(measurement[:,7]+1j*measurement[:,8],dtype="complex")
    S22 = np.array(measurement[:,9]+1j*measurement[:,10],dtype="complex")
    

    Y11= calcul_Y11(S11,S12,S21,S22)
    Y21= calcul_Y21(S11,S12,S21,S22)
    Y12= calcul_Y12(S11,S12,S21,S22)
    Y22= calcul_Y22(S11,S12,S21,S22)
    #! Yij et Zij sont des vecteur contenant les differentes valeurs de Yij/Zij poue different vf, a une freq fixée
    # On fait une moyenne pour retirer la dépendance en freq de Re(Ze)
    Z11,Z12,Z21,Z22=correct_Z_parameter_for_extraction(Y11,Y12,Y21,Y22,freq_indice,open_file,len(vg))
    
    ################################################Extraction partie imaginaire
    
    slope_11 = (Z11[-1].imag  - Z11[0].imag) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
    
    extrapolated_Z11_im = Z11[0].imag + slope_11* (0 - 1/(vg[0]-vth))
    
    slope_12 = (Z12[-1].imag  - Z12[0].imag) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
   
    extrapolated_Z12_im = Z12[0].imag + slope_12 * (0 - 1/(vg[0]-vth))
    
    slope_21= (Z21[-1].imag - Z21[0].imag) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
    
    extrapolated_Z21_im = Z21[0].imag + slope_21 * (0 - 1/(vg[0]-vth))
    
    slope_22 = (Z22[-1].imag - Z22[0].imag) / (1/(vg[-1]-vth) - 1/(vg[0]-vth))
    
    extrapolated_Z22_im = Z22[0].imag + slope_22 * (0 - 1/(vg[0]-vth))
    
    Lg, Ls, Ld = symbols('Xg Xs Xd')
    
    
  
    eq1 = omega*(Ls + Lg) - extrapolated_Z11_im
    eq2 = omega*Ls - extrapolated_Z12_im
    eq3 = omega*(Ld + Ls) - extrapolated_Z22_im

    # Résolution du système d'équations
    solution_L = solve((eq1, eq2, eq3), (Lg, Ls, Ld))
    
    return solution_L[Lg],solution_L[Ls],solution_L[Ld]

#print(extract_Extrinsic_resistance(extract_file,[50,140],open_measure))

#print(extract_Extrinsic_Inductance(extract_file,140,open_measure))





def plot_intrinsic_parameter(file_name,measurement_index,extract_file,open_file): 
    "Plot Z parametre pour un parametre (vg, vbg) fixe selon le fichier"
    "différent d'avant car ici on recupere les valeurs de Parametre pour différentes freq, pas différent Vg"
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
    
    data2=mdm_reader_extract(open_file)
    measurement2=data2[0]
    
    #Open mesure fait pour un seul vg=0 ->on retire le meme open a toute les valeurs de Y11 pour les differents Vg
    S11_op = np.array(measurement2[:,3]+1j*measurement2[:,4],dtype="complex")
    S12_op = np.array(measurement2[:,5]+1j*measurement2[:,6],dtype="complex")
    S21_op = np.array(measurement2[:,7]+1j*measurement2[:,8],dtype="complex")
    S22_op = np.array(measurement2[:,9]+1j*measurement2[:,10],dtype="complex")
    

    Y11_op= calcul_Y11(S11_op,S12_op,S21_op,S22_op)
    Y21_op= calcul_Y21(S11_op,S12_op,S21_op,S22_op)
    Y12_op= calcul_Y12(S11_op,S12_op,S21_op,S22_op)
    Y22_op= calcul_Y22(S11_op,S12_op,S21_op,S22_op)
    
    Y11_tot=Y11-Y11_op
    Y22_tot=Y22-Y22_op
    Y21_tot=Y21-Y21_op
    Y12_tot=Y12-Y12_op
    Y=np.array([[Y11,Y12],[Y21,Y22]])
    
    Rg,Rs,Rd=extract_Extrinsic_resistance(extract_file,[50,140],open_measure)
    Lg,Ls,Ld=extract_Extrinsic_Inductance(extract_file,140,open_measure)
    
    Z_tot=np.zeros((2, 2, len(Y11)), dtype=complex)
    Z_series=np.zeros((2, 2, len(Y11)),dtype=complex)
    for i in range(len(freq)):
        sub_matrix=Y[:,:,i] 
        
        inv = np.linalg.inv(sub_matrix)
        Z_tot[:,:,i]=inv
       
        sub_matrix_2 = np.array([[Rg+Rs+(1j*2*np.pi*freq[i]*(Lg+Ls)), Rs+(1j*2*np.pi*freq[i]*Ls)], [Rs+(1j*2*np.pi*freq[i]*Ls), Rd+Rs+(1j*2*np.pi*freq[i]*(Ld+Ls))]], dtype="complex")
        Z_series[:,:,i]=sub_matrix_2
    
    
        
    Z_device = Z_tot-Z_series
    
    
    
    
    
    
    plt.plot(freq, abs(Z_device[0][0] ), linestyle='-',label="Z11, Vd={},Vg={}".format(vd,vg))
    plt.plot(freq, abs(Z_device[1][0]) , linestyle='-',label="Z21,Vd={},Vg={}".format(vd,vg))
    plt.xlabel('freq[Hz]')
    plt.ylabel('amplitude')
    plt.title('Z parameter after extraction')
    plt.text(1, 20, info[0], fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
    
    plt.legend()
    plt.show()

    
plot_intrinsic_parameter(extract_file,0,extract_file,open_measure)
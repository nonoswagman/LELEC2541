import re
from io import StringIO
import numpy as np

file_name="Die1_SLVTL20_Wf1_Nf16_vbg_0V_idvd_0p8V.mdm"
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

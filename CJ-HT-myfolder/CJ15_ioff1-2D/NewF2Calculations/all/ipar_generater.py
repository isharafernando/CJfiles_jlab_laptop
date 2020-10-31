import os
from os import listdir
from os.path import isfile, join


file_list = [f[:-4] for f in listdir('/users/xjing/QCD_Computation/CJ-code-master/fitpack_offshell/AKP') if f[-4:] == '.par']

print(file_list)

def abbr(filename):
    res = ''
    for i in filename:
        if i not in 'dulea':
            res += i
    return res

for filename in file_list:
    if 'e' in filename:
        new_filename = abbr(filename)
        cmd = 'mv ' + filename + '.par ' + new_filename + '.par'
        print(cmd)
        os.system(cmd)
        cmd = './PDFerr18 ' + new_filename + ' -ipar -deltachi 1.0'
        print(cmd)
        os.system(cmd)



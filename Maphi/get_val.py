import concurrent.futures

#Parsing output file
def check_line(fileout,pattern,num,from_st):
    with open(fileout) as fp:
        for cnt, line in enumerate(fp):
            if cnt > from_st and pattern in line:
                return cnt+int(num)

#here fileout file to parsing, param: search word(s), num: num of line after param, ini: starting line
#fin: end line, from_stt: starying parcing from specific line
def get_val(fileout,param,num,ini,fin,from_stt):
    with open(fileout) as fp:
        for cnt, line in enumerate(fp):
            if cnt == check_line(fileout,param,num,from_stt):
                return line[ini:fin]


with open('test.out') as fp:
    for cnt, line in enumerate(fp):
        if 'SCF CONVERGED AFTER' in line:
            converg=int(cnt)


normal_term = get_val("orca_input.out","****ORCA TERMINATED NORMALLY****",0,0,-1,1)
check_norm="                             ****ORCA TERMINATED NORMALLY****"
conve_term = get_val("orca_input.out","The optimization did not converge",0,0,-1,1)
check_conve="       The optimization did not converge but reached the maximum number of"
if normal_term == check_norm and conve_term == None:
    print("converge")
else:
    print("no converge")
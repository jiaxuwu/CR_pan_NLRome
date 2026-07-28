import sys
codon={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def tr(s):
    s=s.upper().replace('U','T'); out=[]
    for i in range(0,len(s)-2,3):
        out.append(codon.get(s[i:i+3],'X'))
    p=''.join(out)
    if p.endswith('*'): p=p[:-1]
    return p.replace('*','X')  # internal stops -> X so diamond keeps them
inf,outf=sys.argv[1],sys.argv[2]
n=0
with open(inf) as f, open(outf,'w') as o:
    hid=None; seq=[]
    def flush():
        global n
        if hid is not None:
            p=tr(''.join(seq))
            if len(p)>=1:
                o.write('>'+hid+'\n')
                for j in range(0,len(p),60): o.write(p[j:j+60]+'\n')
                n+=1
    for line in f:
        if line.startswith('>'):
            flush(); hid=line[1:].split()[0].strip(); seq=[]
        else: seq.append(line.strip())
    flush()
print(f"{inf}: wrote {n} proteins")

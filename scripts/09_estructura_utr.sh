#!/bin/bash
# Prediccion de estructura secundaria del 5'UTR de la flagelina (RS05045),
# con 100 secuencias control de igual longitud del mismo genoma.
# La flagelina esta en hebra negativa: se toma el complemento reverso.

G=ref/ncbi_dataset/data/GCF_000015445.1/GCF_000015445.1_ASM1544v1_genomic.fna

python3 - "$G" <<'PY' > utr_flagelina.fa
import sys
s=''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>'))
sub=s[1078099:1078240]                       # region intergenica, 141 nt
comp={'A':'T','T':'A','G':'C','C':'G','N':'N'}
print('>5UTR_flagellin_RS05045_141nt')
print(''.join(comp[c] for c in reversed(sub.upper())))
PY

RNAfold --noPS < utr_flagelina.fa

# Control: 100 secuencias aleatorias de igual longitud
python3 - "$G" <<'PY' > control.fa
import sys,random
random.seed(7)
s=''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
for i in range(100):
    p=random.randint(0,len(s)-142)
    print(f'>ctrl{i}'); print(s[p:p+141])
PY

RNAfold --noPS < control.fa | grep -oP '\(\s*-?[0-9]+\.[0-9]+\)$' | tr -d '() ' > mfe_ctrl.txt
python3 -c "
import statistics as st
v=sorted(float(x) for x in open('mfe_ctrl.txt'))
print('Control: media %.1f  mediana %.1f  percentil5 %.1f'%(st.mean(v),st.median(v),v[5]))"

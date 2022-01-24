import re

lines = []
with open('data/raw_snps.vcf', 'r') as v:
    for l in v:
        lines.append(l)

und = []
ref = []
het = []
alt = []

for v in lines:
    und.append(re.findall('\\./\\.', v))
    ref.append(re.findall('0/0', v))
    het.append(re.findall('0/1', v))
    alt.append(re.findall('1/1', v))

lund = [len(x) for x in und]
lref = [len(x) for x in ref]
lhet = [len(x) for x in het]
lalt = [len(x) for x in alt]

counts = list(map(lambda x, y, z, t: x + y + z + t, lund, lref, lhet, lalt))

vcf = []
for i, j in enumerate(lines):
    if i in list(range(39)) or counts[i] == 159:
        vcf.append(j)

output = open('data/biallelic_snps.vcf', 'w')
for i in range(len(vcf)):
    output.write(vcf[i])
output.close()

from gwf import Workflow

# Templates -------------------------------------
def quality_filter(inb, oub):
    inputs = [inb]
    outputs = [oub]
    options = {
        'cores': 2,
        'memory': '128g',
        'walltime': '12:00:00'
    }
    spec = '''
samtools view -bq 30 {inb} > {oub}
    '''.format(inb=inb, oub=oub)
    return inputs, outputs, options, spec


def index(bam, bai):
    inputs = [bam]
    outputs = [bai]
    options = {
        'cores': 2,
        'memory': '128g',
        'walltime': '12:00:00'
    }
    spec = '''
samtools index {bam} {bai}
    '''.format(bam=bam, bai=bai)
    return inputs, outputs, options, spec


# Workflow --------------------------------------
files = []
with open('bam_files.txt', 'rt') as f:
    for line in f:
        files.append(line.strip())
accessions = [x[5:22] for x in files]

gwf = Workflow(defaults={'account': 'InRoot'})

for i, j in enumerate(files):
    gwf.target_from_template('quality_{}'.format(accessions[i]),
    quality_filter(
        inb=j,
        oub='filtered_data/{}.bam'.format(accessions[i])
    ))
    
    gwf.target_from_template('index_{}'.format(accessions[i]),
                            index(
                                     bam='filtered_data/{}.bam'.format(accessions[i]),
                                     bai='filtered_data/{}.bam.bai'.format(accessions[i])
                                 ))

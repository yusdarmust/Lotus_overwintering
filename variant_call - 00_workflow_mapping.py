from gwf import Workflow
import glob

# Templates -------------------------------------
def bwa_map(fa, fq1, fq2, output):
    '''
    Template for mapping short reads to the genome.
    '''
    inputs = [fa, fq1, fq2,
              '{}.amb'.format(fa),
              '{}.ann'.format(fa),
              '{}.bwt'.format(fa),
              '{}.pac'.format(fa),
              '{}.sa'.format(fa)]
    outputs = [output]
    options = {
        'cores': 24,
        'memory': '128g',
        'walltime': '240:00:00'
    }
    spec = '''
bwa mem {fa} {fq1} {fq2} | samtools view -Sb > {output}
'''.format(fa=fa, fq1=fq1, fq2=fq2, output=output)
    return inputs, outputs, options, spec


def samtools_sort(mapped, sorted_):
    '''
    Template for sorting the mapped reads.
    '''
    inputs = [mapped]
    outputs = [sorted_]
    options = {
        'cores': 12,
        'memory': '64g',
        'walltime': '120:00:00'
    }
    spec = '''
samtools sort {mapped} > {sorted_}
    '''.format(mapped=mapped, sorted_=sorted_)
    return inputs, outputs, options, spec


def samtools_index(bam, bai):
    '''
    Template for indexing the sorted reads.
    '''
    inputs = [bam]
    outputs = [bai]
    options = {
        'cores': 12,
        'memory': '64g',
        'walltime': '120:00:00'
    }
    spec = '''
samtools index {bam} {bai}
    '''.format(bam=bam, bai=bai)
    return inputs, outputs, options, spec

# Workflow --------------------------------------
files1 = sorted(glob.glob('data/*_R1.fastq'))
files2 = sorted(glob.glob('data/*_R2.fastq'))
files = list(zip(files1, files2))
rg = '../../InRoot/Backup/data/02_20200402_Gifu1.2_ref_data/LjGifu1.1_pseudomol.fa'

gwf = Workflow()

for i in range(len(files)):
        gwf.target_from_template('bwaMapping_{}'.format(files[i][0][5:10]),
                                 bwa_map(fa=rg,
                                         fq1=files[i][0],
                                         fq2=files[i][1],
                                         output='results/mapped_{}.bam'.format(
                                             files[i][0][5:10])
                                         ))

        gwf.target_from_template('samtoolsSort_{}'.format(files[i][0][5:10]),
                                 samtools_sort(
                                     mapped='results/mapped_{}.bam'.format(
                                         files[i][0][5:10]),
                                     sorted_='results/sorted_{}.bam'.format(
                                         files[i][0][5:10])
                                     ))

        gwf.target_from_template('samtoolsIndex_{}'.format(files[i][0][5:10]),
                                 samtools_index(
                                     bam='results/sorted_{}.bam'.format(
                                         files[i][0][5:10]),
                                     bai='results/sorted_{}.bam.bai'.format(
                                         files[i][0][5:10])
                                     ))

from gwf import Workflow

# Templates -------------------------------------
def variant_calling(fa, bamlist, output):
    inputs = [fa, bamlist]
    outputs = [output]
    options = {
        'cores': 16,
        'memory': '128g',
        'walltime': '160:00:00'
    }
    spec = '''
bcftools mpileup -Ou -f {fa} -b {bamlist} |
bcftools call -mv > {output}
    '''.format(fa=fa, bamlist=bamlist, output=output)
    return inputs, outputs, options, spec

# Workflow --------------------------------------
chromosome = ['LjG1.1_chr1', 'LjG1.1_chr2', 'LjG1.1_chr3',
              'LjG1.1_chr4', 'LjG1.1_chr5', 'LjG1.1_chr6',
              'LjG1.1_chr0', 'LjG1.1_chloro', 'LjG1.1_mito']    
rg = '../../InRoot/Backup/data/02_20200402_Gifu1.2_ref_data/LjGifu1.1_pseudomol.fa'

gwf = Workflow(defaults={'account': 'InRoot'})

for c in chromosome:
    gwf.target_from_template('variant_calling_{}'.format(c),
                             variant_calling(
                                 fa=rg,
                                 bamlist='filtered_data/' + c + '.txt',
                                 output='filtered_results/variant_calling_{}.vcf'.format(c)
                             ))

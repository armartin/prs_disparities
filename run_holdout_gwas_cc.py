from __future__ import print_function
import sys
# from gnomad_hail.utils.slack_utils import *
from hail import *
from pprint import pprint
import argparse
hc = HailContext()

def main(args):
    print('Read/assemble phenotypes')
    """
    need to get these from a bunch of split up tables
    if it's not possible to get these from hail easily, potentially read them in with pandas,
    where you can scan the header before deciding which phenotypes to read
    """
    sex = 'both_sexes'
    kt_phenotype = hc.import_table('gs://mkanai/disparities/ukb31063.phecode_5diseases.both_sexes.tsv.bgz',
                                    key='s', impute=True, types={'s': TString()}, missing='NA')
    phenos = ['crc', 't2d', 'glaucoma', 'afib', 'ra']
    renamed = {'s': 's', 'CRC': 'crc', 'T2D': 't2d', 'Glaucoma': 'glaucoma', 'AFib': 'afib', 'RA': 'ra'}
    kt_phenotype = kt_phenotype.rename(renamed)

    pheno_gwas_cases = {'crc': 3461, 't2d': 15303, 'glaucoma': 3383, 'afib': 7674, 'ra': 3297}
    pheno_gwas_controls = {'crc': 153621, 't2d': 130904, 'glaucoma': 164288, 'afib': 153907, 'ra': 165386}
    pheno_holdout_cases = {'crc': 500, 't2d': 500, 'glaucoma': 500, 'afib': 500, 'ra': 500}
    pheno_holdout_controls = {'crc': 500, 't2d': 500, 'glaucoma': 500, 'afib': 500, 'ra': 500}

    print('Read covariates, QC sites')
    vds_variants = hc.read('gs://phenotype_31063/hail/0.1/ukb31063.gwas_variants.autosomes.with_qc_annotations.vds')
    kt_covariates = hc.read_table('gs://phenotype_31063/hail/0.1//ukb31063.gwas_covariates.{}.kt'.format(sex))
    # 655633986626-compute@developer.gserviceaccount.com does not have storage.buckets.get access to phenotype_31063.

    print('Reading UKBB imputed variants')
    import_expr = '{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
    # import_expr = '{22}'

    covariate_expr = ['sa.covariates.age',
                      'sa.covariates.age_squared',
                      'sa.covariates.isFemale',
                      'sa.covariates.age_isFemale',
                      'sa.covariates.age_squared_isFemale'] + \
                     ['sa.covariates.PC{:}'.format(i) for i in xrange(1, 21)]

    vds = (
    hc.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'.format(import_expr),
                   sample_file='gs://phenotype_31063/ukb31063.{}.sample'.format('autosomes'),
                   tolerance=0.2)
        .annotate_variants_vds(vds_variants, expr='va.AF = vds.qc.AF, va.info = vds.info')
        .filter_variants_expr('isDefined(va.AF)', keep=True)
        .annotate_samples_table(kt_covariates, root='sa.covariates')
        .annotate_samples_table(kt_phenotype, root='sa.phenotypes')
    )

    vds.variants_table().write('gs://mkanai/disparities/pheno_31063_holdout_gwas_cc.kt', args.overwrite)
    kt_results = hc.read_table('gs://mkanai/disparities/pheno_31063_holdout_gwas_cc.kt')

    print('Run regressions per phenotype')
    for i in range(3, len(phenos)):
        pheno_code = phenos[i]
        #print(pheno_codes[i])
        """
        0. Subset to samples whose phenotypes are not missing
        1. Add a random column to my key table: kt_phenotype.annotate('rand=rnorm(0,1)')
        2. Order by this random column: kt_phenotype.order_by('rand')
        3. Add the index: kt_phenotype.indexed()
        4. Keep GWAS inds from the smallest rand, holdout inds from the next smallest, forget everyone else
        """

        pheno_subset = vds.filter_samples_expr('isDefined(sa.phenotypes.`{}`) && '
                                               'isDefined(sa.covariates.PC1)'.format(pheno_code))
        pheno_controls = vds.filter_samples_expr('sa.phenotypes.`{}` == 0'.format(pheno_code), keep=True).num_samples
        print(pheno_subset.count())
        pheno_sub_table = pheno_subset.samples_table()
        pheno_sub_table = pheno_sub_table.annotate('rand=10*sa.phenotypes.`{}`+runif(0,1)'.format(pheno_code))
        pheno_sub_table = pheno_sub_table.order_by('rand')
        pheno_sub_table = pheno_sub_table.indexed()
        pheno_sub_table = pheno_sub_table.annotate('gwas_holdout= if (index < {}) "gwas" '
                                                   'else if (index < {}) "holdout" '
                                                   'else if (index < {} && rand >= 10) "gwas" '
                                                   'else if (index < {} && rand >= 10) "holdout" '
                                                   'else "NA"'.format(pheno_gwas_controls[pheno_code],
                                                                      pheno_gwas_controls[pheno_code] +
                                                                      pheno_holdout_controls[pheno_code],
                                                                      pheno_controls +
                                                                      pheno_gwas_cases[pheno_code],
                                                                      pheno_controls +
                                                                      pheno_gwas_cases[pheno_code] +
                                                                      pheno_holdout_cases[pheno_code]))
        pheno_sub_print = pheno_sub_table.drop(['sa', 'rand', 'index'])
        pheno_sub_print.export('gs://mkanai/disparities/pheno_31063_holdout_gwas_{}.info.txt.gz'.format(pheno_code))
        pheno_sub_print = hc.import_table('gs://mkanai/disparities/pheno_31063_holdout_gwas_{}.info.txt.gz'.format(pheno_code),
                                          key='s', impute=True, types={'s': TString()}, missing='NA')

        my_vds = vds.filter_samples_table(pheno_sub_table.filter('gwas_holdout=="gwas"'))

        my_vds = my_vds.logreg('wald', y='sa.phenotypes.`{}`'.format(pheno_code),
                             covariates=covariate_expr,
                             use_dosages=True)
        # my_vds = my_vds.annotate_variants_expr('va.results' = ['va.linreg'])
        kt_results = my_vds.variants_table()
        kt_export = kt_results.annotate(['chr = v.contig',
                                         'pos = v.start',
                                         'ref = v.ref',
                                         'alt = v.alt',
                                         'rsid = va.rsid',
                                      #  'nCompleteSamples = va.linreg.nCompleteSamples',
                                      #  'AC = va.linreg.AC',
                                      #  'ytx = va.linreg.ytx[0]',
                                         'beta = va.logreg.beta',
                                         'se = va.logreg.se',
                                         'zstat = va.logreg.zstat',
                                         'pval = va.logreg.pval'])

        kt_export2 = kt_export.drop(['v', 'va'])
        kt_export2.write('gs://mkanai/disparities/pheno_31063_holdout_gwas_{}.kt'.format(pheno_code), args.overwrite)
        kt_export2 = hc.read_table('gs://mkanai/disparities/pheno_31063_holdout_gwas_{}.kt'.format(pheno_code))
        kt_export2.export('gs://mkanai/disparities/pheno_31063_holdout_gwas_{}.txt.gz'.format(pheno_code))


if __name__ == '__main__':
    print('Starting run')
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
    #  try_slack('@armartin', main, args)

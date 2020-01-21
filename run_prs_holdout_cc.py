import hail as hl
import pickle
import time
import argparse
from pprint import pprint

def flip_text(base):
    """
    :param StringExpression base: Expression of a single base
    :return: StringExpression of flipped base
    :rtype: StringExpression
    """
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('C', 'G')
            .when('G', 'C')
            .default(base))


def annotate_beta(mt, ss_loc):
    mt = mt.annotate_rows(**{
        'beta': hl.case()
                          .when(((mt.alleles[0] == ss_loc.ref) &
                                 (mt.alleles[1] == ss_loc.alt)) |
                                ((flip_text(mt.alleles[0]) == ss_loc.ref) &
                                 (flip_text(mt.alleles[1]) == ss_loc.alt)),
                                ss_loc.beta)
                          .when(((mt.alleles[0] == ss_loc.alt) &
                                 (mt.alleles[1] == ss_loc.ref)) |
                                ((flip_text(mt.alleles[0]) == ss_loc.alt) &
                                 (flip_text(mt.alleles[1]) == ss_loc.ref)),
                                (-1 * ss_loc.beta))
                          .or_missing()}
                          )
    return(mt)


def specific_clumps(filename):
    clump = hl.import_table(filename, delimiter='\s+', min_partitions=10, types={'P': hl.tfloat}, skip_blank_lines=True)
    clump = clump.key_by(locus = hl.locus(hl.str(clump.CHR), hl.int(clump.BP)))
    return clump


def main(args):
    ########################################################################
    ### initialize
    phenos = ['crc', 't2d', 'glaucoma', 'afib', 'ra']
    renamed = {'s': 's', 'CRC': 'crc', 'T2D': 't2d', 'Glaucoma': 'glaucoma', 'AFib': 'afib', 'RA': 'ra'}
    phenotype = 'ALL5cc'
    sumstats_text_file = args.dirname + args.basename + 'ALL5cc.clumped'
    prs_loci_table_location = args.dirname + 'keytables/ukb-'+phenotype+'-pt-sumstats-locus-allele-keyed.kt'
    contig_row_dict_location = args.dirname + 'contig_row_dict-'+phenotype

    contigs = {'0{}'.format(x):str(x) for x in range(1, 10)}

    bgen_files = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}_v3.bgen'

    start = time.time()
    # large block size because we read very little data (due to filtering & ignoring genotypes)
    # hl.init(branching_factor=10, min_block_size=2000)
    hl.init()


    ################################################################################
    ### set up the sumstats table (chr, bp for union SNPs)
    if (args.generate_prs_loci_table):
        t = hl.import_table(sumstats_text_file,
                            delimiter='\s+',
                            impute=True)
        t = t.select(locus = hl.locus(hl.str(t.CHR), t.BP))
        t = t.key_by('locus')
        t.write(prs_loci_table_location, overwrite=True)

    ss = hl.read_table(prs_loci_table_location)

    ################################################################################
    ### Get true phenotypes from UKBB
    if args.pheno_table:
        phenotypes = hl.import_table('gs://mkanai/disparities/ukb31063.phecode_5diseases.both_sexes.tsv.bgz',
                                     key='s', impute=True, types={'s': hl.tstr})
        phenotypes = phenotypes.rename(renamed)

        covariates = hl.import_table('gs://phenotype_31063/ukb31063.gwas_covariates.both_sexes.tsv',
                                     key='s', impute=True, types={'s': hl.tstr})

        samples = covariates.annotate(**phenotypes[covariates.s])

        # Write pheno/covar/sample info table
        for pheno in phenos:
            gwas_holdout = hl.import_table('gs://mkanai/disparities/ukbb/pheno_31063_holdout_gwas_' + pheno + '.info.txt.gz', delimiter='\s+').key_by('s')

            samples = samples.annotate(**{pheno + '_holdout': gwas_holdout[samples.s].gwas_holdout == 'holdout'})

        samples.write('gs://mkanai/disparities/pheno_31063_holdout_gwas_cc_phenos.ht', args.overwrite)

    if args.ss_tables:
        # Write ss info
        for pheno in phenos:
            print(pheno)
            ss = hl.import_table(args.dirname + args.basename + pheno + '.*.bgz',
                                 delimiter='\s+',
                                 impute=True,
                                 types={'beta': hl.tfloat, 'pval': hl.tfloat, 'pos': hl.tint,
                                        'nCompleteSamples': hl.tint, 'AC': hl.tfloat, 'ytx': hl.tfloat, 'se': hl.tfloat,
                                        'tstat': hl.tfloat})
            ss = ss.key_by(locus = hl.locus(hl.str(ss.chr), hl.int(ss.pos))).repartition(200)

            ss.write(args.dirname + args.basename + pheno + '.ht', True)

    ################################################################################
    ### Run the PRS using phenotype-specific clump variants
    if args.write_bgen:
        mt_all = hl.import_bgen(
            bgen_files,
            entry_fields=['dosage'],
            sample_file='gs://phenotype_31063/ukb31063.autosomes.sample',
            variants=ss.locus)

        samples = hl.read_table('gs://mkanai/disparities/pheno_31063_holdout_gwas_cc_phenos.ht')
        mt_all = mt_all.annotate_cols(**samples[mt_all.s]) # ok that phenos keyed on userId not s?

        mt_all.repartition(5000, shuffle=False).write(args.dirname + args.basename + 'ALL5cc.mt', args.overwrite)

    mt_all = hl.read_matrix_table(args.dirname + args.basename + 'ALL5cc.mt')


    for pheno in phenos: #[6:len(phenos)]:
        print(pheno)
        ss = hl.read_table(args.dirname + args.basename + pheno + '.ht')

        """
        To add:
        - Filter only to samples in holdout GWAS
        - Filter to rows in phenotype-specific clump file
        - Build PRS for 10 p-value thresholds
        - Also fix nt1/nt2 to A1 and A2 (check) from sumstats.
        """
        # filter to only samples held out from GWAS
        mt = mt_all.filter_cols(mt_all[pheno + '_holdout'])

        mt = mt.annotate_rows(ss=ss[mt.locus])
        mt = annotate_beta(mt, mt.ss)

        # p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5, 's10': 1}
        p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2}

        pheno_clump = specific_clumps(args.dirname + args.basename + pheno + '.clumped')

        mt = mt.filter_rows(hl.is_defined(pheno_clump[mt.locus]))
        print(mt.count())

        annot_expr = {
            k: hl.agg.sum(mt.beta * mt.dosage * hl.int(mt.ss.pval < v))
            for k, v in p_max.items()}

        mt = mt.annotate_cols(**annot_expr)

        mt.cols().write(args.dirname + 'UKB_' + pheno + '_PRS.ht', stage_locally=True, overwrite=True)
        ht = hl.read_table(args.dirname + 'UKB_' + pheno + '_PRS.ht')
        ht_out = ht.drop(*[x for x in list(ht.row) if 'holdout' in x], *[x for x in phenos if pheno not in x])

        output_location = args.dirname + 'UKB_' + pheno + '_PRS.txt.bgz'
        ht_out.export(output_location)
    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--generate_prs_loci_table', action='store_true')
    parser.add_argument('--pheno_table', action='store_true')
    parser.add_argument('--ss_tables', action='store_true')
    parser.add_argument('--write_bgen', action='store_true')
    parser.add_argument('--dirname', default='gs://mkanai/disparities/ukbb/') # gs://armartin/disparities/ukbb/
    parser.add_argument('--basename', default='pheno_31063_holdout_gwas_')  # pheno_31063_holdout_gwas_

    args = parser.parse_args()
    main(args)
    #try_slack('@armartin', main, args)

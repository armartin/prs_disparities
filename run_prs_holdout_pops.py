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

def specific_clumps(filename):
    clump = hl.import_table(filename, delimiter='\s+', min_partitions=10, types={'P': hl.tfloat})
    clump_dict = clump.aggregate(hl.dict(hl.agg.collect(
        (hl.locus(hl.str(clump.CHR), hl.int(clump.BP)),
        True)
    )))
    return hl.literal(clump_dict)

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


def main(args):
    ########################################################################
    ### initialize
    phenos = ['height', 'bmi', 'sbp', 'dbp', 'wbc', 'monocyte', 'neutrophil', 'eosinophil', 'basophil', 'lymphocyte',
              'rbc', 'mch', 'mcv', 'mchc', 'hb', 'ht', 'plt']
    phenotype = 'ALL17'
    sumstats_text_file = args.dirname + args.basename + 'ALL17.clumped'
    prs_loci_table_location = args.dirname + 'keytables/ukb-'+phenotype+'-pt-sumstats-locus-allele-keyed.kt'
    contig_row_dict_location = args.dirname + 'contig_row_dict-'+phenotype

    contigs = {'0{}'.format(x):str(x) for x in range(1, 10)}

    bgen_files = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}_v3.bgen'

    start = time.time()
    # large block size because we read very little data (due to filtering & ignoring genotypes)
    hl.init(branching_factor=10, min_block_size=2000)

    mt_all = hl.read_matrix_table(args.dirname + args.basename + 'ALL17.mt')


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
        
        mt = mt_all.annotate_rows(ss=ss[mt_all.locus])
        mt = annotate_beta(mt, mt.ss)

        p_max = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5, 's10': 1}

        pheno_clump = specific_clumps(args.dirname + args.basename + pheno + '.clumped')

        mt = mt.filter_rows(pheno_clump.get(mt.locus, False))
        print(mt.count())

        annot_expr = {
            k: hl.agg.sum(mt.beta * mt.dosage * hl.int(mt.ss.pval < v))
            for k, v in p_max.items()}

        mt = mt.annotate_cols(**annot_expr)

        mt.cols().write(args.outdir + 'UKB_' + pheno + '_PRS.ht', stage_locally=True, overwrite=True)
        ht = hl.read_table(args.outdir + 'UKB_' + pheno + '_PRS.ht')
        ht_out = ht.drop(*[x for x in list(ht.row) if 'holdout' in x], *[x for x in phenos if pheno not in x])

        output_location = args.outdir + 'UKB_' + pheno + '_PRS.txt.bgz'
        ht_out.export(output_location)
    end = time.time()
    print("Success! Job was completed in %s" % time.strftime("%H:%M:%S", time.gmtime(end - start)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--generate_prs_loci_table', action='store_true')
    parser.add_argument('--generate_contig_row_dict', action='store_true')
    parser.add_argument('--pheno_table', action='store_true')
    parser.add_argument('--ss_tables', action='store_true')
    parser.add_argument('--write_bgen', action='store_true')
    parser.add_argument('--dirname', default='gs://armartin/disparities/bbj/') # gs://armartin/disparities/ukbb/
    parser.add_argument('--outdir', default='gs://armartin/disparities/bbj/')  # gs://armartin/disparities/ukbb/
    parser.add_argument('--basename', default='BBJ_holdout_gwas_')  # pheno_31063_holdout_gwas_

    args = parser.parse_args()
    main(args)
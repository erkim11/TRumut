/*
main.cpp: `trcall` subcommand dispatcher.
*/

#include <stdio.h>
#include <string.h>

#ifndef TRCALL_VERSION
#define TRCALL_VERSION "0.1.0"
#endif

int allele_main(int argc, char *argv[]);
int prep_catalog_main(int argc, char *argv[]);
int build_pon_main(int argc, char *argv[]);
int use_pon_main(int argc, char *argv[]);
int cancer_main(int argc, char *argv[]);
int trio_main(int argc, char *argv[]);
int mutual_main(int argc, char *argv[]);

static void usage(FILE *fp)
{
    fprintf(fp,
"trcall %s -- tandem repeat allele calling\n"
"\n"
"Usage: trcall <command> [options]\n"
"\n"
"Main workflow:\n"
"  allele        call tandem repeat alleles per sample; input to the three below\n"
"  cancer        alleles in a tumor that are absent from its matched normal\n"
"  trio          alleles in a child that are absent from both parents\n"
"  use-pon       alleles outside the range seen in a panel of normals, for samples with no matched control\n"
"\n"
"Not needed for the main workflow:\n"
"  prep-catalog  build a uniqmer catalog for allele from your own tandem repeat catalog (a pre-built one exists)\n"
"  build-pon     build a panel of normals for use-pon from your own samples (a pre-built one exists)\n"
"  mutual        compare two or more samples symmetrically: alleles present in one but absent from another\n"
"\n"
"  trcall <command> -h   options of one command\n"
"  -h, --help            show this message\n"
"  -v, --version         show the version\n"
"\n"
"Documentation: https://erkim11.github.io/trcall/\n",
            TRCALL_VERSION);
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(stderr); return 1; }

    const char *cmd = argv[1];

    if (strcmp(cmd, "-h") == 0 || strcmp(cmd, "--help") == 0 || strcmp(cmd, "help") == 0) {
        usage(stdout);
        return 0;
    }
    if (strcmp(cmd, "-v") == 0 || strcmp(cmd, "--version") == 0) {
        printf("trcall %s\n", TRCALL_VERSION);
        return 0;
    }

    char progname[64];
    snprintf(progname, sizeof progname, "trcall %s", cmd);
    char *saved = argv[1];
    argv[1] = progname;

    int rc;
    if (strcmp(cmd, "allele") == 0) {
        rc = allele_main(argc - 1, argv + 1);
    } else if (strcmp(cmd, "prep-catalog") == 0) {
        rc = prep_catalog_main(argc - 1, argv + 1);
    } else if (strcmp(cmd, "build-pon") == 0) {
        rc = build_pon_main(argc - 1, argv + 1);
    } else if (strcmp(cmd, "use-pon") == 0) {
        rc = use_pon_main(argc - 1, argv + 1);
    } else if (strcmp(cmd, "cancer") == 0) {
        rc = cancer_main(argc - 1, argv + 1);
    } else if (strcmp(cmd, "trio") == 0) {
        rc = trio_main(argc - 1, argv + 1);
    } else if (strcmp(cmd, "mutual") == 0) {
        rc = mutual_main(argc - 1, argv + 1);
    } else {
        argv[1] = saved;
        fprintf(stderr, "[ERROR] unknown command '%s'\n\n", cmd);
        usage(stderr);
        return 1;
    }

    argv[1] = saved;
    return rc;
}

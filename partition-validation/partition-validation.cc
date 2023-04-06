#include "partition-validation.h"

using namespace std;

//unsigned stou(char *s) {
//    return strtoul(s, (char **) NULL, 10);
//}

// Call: trade <seed> <Ntries>
int main(int argc, char *argv[]) {

    cout << "Version: Nov 5, 2022.\n";
    cout << "Command: ";
    cout << argv[0];
    for (int i = 1; i < argc; i++)
        cout << " " << argv[i];
    cout << '\n';

    // Parse command input
    const string CALL_SYNTAX = "Call: ./partition-validation [-h] [-s <seed>] [-t <distance threshold>] "
                               "[--validate <validation size>] [--k-fold-crossvalidate k] "
                               "[--validation-sampling <training size> <validation size> <validation samples>] "
                               "input_partitions.txt output_clustering_txt\n";
    if (argc == 1) {
        cout << CALL_SYNTAX;
        exit(-1);
    }

    int seed = 1234;

    string inFileName = "noname";
    string outFileName = "noname";

    int argNr = 1;
    double distThreshold = 0.2;
    int nValidationPartitions = 0;
    int nTrainingPartitions = 0;
    int crossValidateK = 0;
    int validationSamples = 0;

    while (argNr < argc) {
        if (to_string(argv[argNr]) == "-h") {
            cout << CALL_SYNTAX;
            cout << "-h: This help\n";
            cout << "seed: Any positive integer.\n";
            cout << "distance threshold: "
                    "The max distance from center partition to other partitions in the cluster. Default: 0.2.\n";
            cout << "--validate <validation size>: "
                    "The number of partitions at the end of the input that will be used for validation. "
                    "The first partitions will be used to find clusters. Default is 0 validation partitions.\n";
            cout << "--k-fold-crossvalidate k: "
                    "Perform k-fold cross-validation of all partitions. The training partitions will be used to "
                    "find clusters and the other ones for validation. Default is 0 folds for no cross-validation.\n";
            cout << "--validation-sampling <training size> <validation size> <validation samples>: "
                    "A random set of validation size partitions will be held out from training size partitions "
                    "in each of validation samples resamplings. Reports the average fraction of validation partitions "
                    "that belong to clusters. Default is no validation sampling.\n";
            cout << "input_partitions.txt: "
                    "Each column corresponds to a partition and each row corresponds to a node id.\n";
            cout << "output_clustering.txt: clusterID partitionID\n";
            exit(-1);
        } else if (to_string(argv[argNr]) == "-s") {
            argNr++;
            seed = atoi(argv[argNr]);
            argNr++;
        } else if (to_string(argv[argNr]) == "-t") {
            argNr++;
            distThreshold = atof(argv[argNr]);
            argNr++;
        } else if (to_string(argv[argNr]) == "--validate") {
            argNr++;
            nValidationPartitions = atoi(argv[argNr]);
            argNr++;
        } else if (to_string(argv[argNr]) == "--k-fold-crossvalidate") {
            argNr++;
            crossValidateK = atoi(argv[argNr]);
            argNr++;
        } else if (to_string(argv[argNr]) == "--validation-sampling") {
            argNr++;
            nTrainingPartitions = atoi(argv[argNr]);
            argNr++;
            nValidationPartitions = atoi(argv[argNr]);
            argNr++;
            validationSamples = atoi(argv[argNr]);
            argNr++;
        } else {
            if (argv[argNr][0] == '-') {
                cout << "Unknown command: " << to_string(argv[argNr]) << '\n';
                cout << CALL_SYNTAX;
                exit(-1);
            }

            if (inFileName == "noname") {
                inFileName = string(argv[argNr]);
                argNr++;
            } else if (outFileName == "noname") {
                outFileName = string(argv[argNr]);
                argNr++;
            }
        }
    }

    if (inFileName == "noname") {
        cout << "Missing infile\n";
        cout << CALL_SYNTAX;
        exit(-1);
    }
    if (outFileName == "noname") {
        cout << "Missing outfile\n";
        cout << CALL_SYNTAX;
        exit(-1);
    }

    cout << "Setup:\n";
    cout << "--> Using seed: " << seed << '\n';
    cout << "--> Will cluster partitions such that no partition is farther away from its center than: " << distThreshold
         << '\n';
    cout << "--> Will read partitions from file: " << inFileName << '\n';

    if (crossValidateK > 0) {
        cout << "--> Will perform " << crossValidateK << "-fold cross-valiation.\n";
        nValidationPartitions = 0;
        validationSamples = 0;
    } else if (validationSamples > 0) {
        cout << "--> Will perform validation sampling of " << nValidationPartitions << " partitions given "
             << nTrainingPartitions << " training partitions " << validationSamples << " times.\n";
        crossValidateK = 0;
    } else if (nValidationPartitions > 0) {
        cout << "--> Will use the last " << nValidationPartitions << " partitions for validation.\n";
        crossValidateK = 0;
        validationSamples = 0;
    }

    cout << "--> Will write clusters to file: " << outFileName << '\n';
    cout << "--> Will use number of threads: " << omp_get_max_threads() << '\n';

    Partitions partitions(inFileName, outFileName, distThreshold, nTrainingPartitions, nValidationPartitions,
                          crossValidateK, validationSamples, seed);

    if (crossValidateK == 0) {
        if (validationSamples == 0) {
            partitions.clusterPartitions(0);
            if (nValidationPartitions > 0) {
                partitions.validatePartitions(0);
                cout << "--> Fraction of validation partitions that fits in a cluster: "
                     << 1.0 * partitions.nTotValidated / partitions.nTotTested << '\n';
            }
            partitions.printClusters();
        } else {
            for (int sample = 0; sample < validationSamples; sample++) {
                partitions.clusterPartitions(0);
                partitions.validatePartitions(0);
                cout << "--> Fraction of validation partitions that fits in a cluster after " << sample + 1
                     << " samples: " << 1.0 * partitions.nTotValidated / partitions.nTotTested << '\n';
            }
        }
    } else {
        for (int fold = 0; fold < crossValidateK; fold++) {
            cout << '\n' << "Fold " << fold + 1 << "/" << crossValidateK << '\n';
            partitions.clusterPartitions(fold);
            partitions.validatePartitions(fold);
            cout << "--> Fraction of validation partitions that fits in a cluster after " << fold + 1 << " folds: "
                 << 1.0 * partitions.nTotValidated / partitions.nTotTested << '\n';
        }
    }
}
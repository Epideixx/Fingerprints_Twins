import argparse
from PSD_accuracy import main as compute_accuracy
from PSD_accuracy_30_sec import main as compute_accuracy_30_sec
from PSD_accuracy_empty_room import main as compute_accuracy_empty_room
from PSD_correlations import main as compute_correlations
from PSD_Heritability import main as compute_heritability
from PSD_Heritability_Anatomy import main as compute_heritability_anatomy
from PSD_ICC_Fingerprint import main as compute_icc_fingerprint


parser = argparse.ArgumentParser(description='PyTorch PennTreeBank RNN/LSTM Language Model')
parser.add_argument('--acc', type=bool, default=True,
                    help='either if we compute the accuracy or not')
parser.add_argument('--acc_30_sec', type=bool, default=False,
                    help='either if we compute the accuracy or not for the 30-sec PSDs')
parser.add_argument('--acc_empty_room', type=bool, default=False,
                    help='either if we compute the accuracy of the empty room or not')
parser.add_argument('--corr', type=bool, default=False,
                    help='either if we compute the correlations or not')
parser.add_argument('--heritability', type=bool, default=True,
                    help='either if we compute the heritability or not')
parser.add_argument('--heritability_anatomy', type=bool, default=False,
                    help='either if we compute the hertiability of anatomy or not')
parser.add_argument('--fingerprint', type=bool, default=True,
                    help='either if we compute the ICC for fingerprinting or not')
parser.add_argument('--data', type=str, default='Data/Schaefer',
                    help='location of the data')
parser.add_argument('--data_30_sec', type=str, default='Data/Schaefer_30_second',
                    help='location of the data for 30-seconds PSDs')
parser.add_argument('--results', type=str, default='Results_Log_Schaefer_test',
                    help='location of the results')
parser.add_argument('--only_gt', type=bool, default=True,
                    help='choice of if we use only twin pairs based on genetic test (True), or if we also use the self reported twins (False)')
parser.add_argument('--correlation_type', type=str, default="icc",
                    help='correlation used for heritability (pearson or icc)')
parser.add_argument('--n_resample', type=int, default=1000,
                    help='number of bootstraps')
parser.add_argument('--icc_without_twins', type=bool, default=False,
                    help='either ICC for fingerprinting is computed with or without twins')
args = parser.parse_args()

print(args)

if args.acc:
    compute_accuracy(data_path=args.data, main_folder_results=args.results, only_gt=args.only_gt, n_resample=args.n_resample)

if args.acc_30_sec:
    compute_accuracy_30_sec(data_path=args.data_30_sec, main_folder_results=args.results, only_gt=args.only_gt)

if args.acc_empty_room:
    compute_accuracy_empty_room(data_path=args.data, main_folder_results=args.results, only_gt=args.only_gt, n_resample=args.n_resample)

if args.corr:
    compute_correlations(data_path=args.data, main_folder_results=args.results, only_gt=args.only_gt)

if args.heritabiliy:
    compute_heritability(data_path=args.data, main_folder_results=args.results, only_gt=args.only_gt, correlation_type = args.correlation_type)

if args.heritability_anatomy :
    compute_heritability_anatomy(data_path=args.data, main_folder_results=args.results, only_gt=args.only_gt, correlation_type = args.correlation_type)

if args.fingerprint:
    compute_icc_fingerprint(data_path=args.data, main_folder_results=args.results, only_gt=args.only_gt, n_resample=args.n_resample, without_twins=args.icc_without_twins)
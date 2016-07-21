import os
import math

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
        io_path = '../mixed_effects_results/'
        hypothesis_test_path = '/home/sci/prasanna/work/git/PHD_ShapeGrant/shape-lme/bin/hypothesis-testing-tutorial'
        template_file_dir = '/home/sci/prasanna/work/git/PHD_ShapeGrant/deformetrica/tutorial/shape_lpts_for_stats/50855/'

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	times = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]
        hd = ['3', '3', '3', '0', '0', '3', '0', '3', '3', '3', '0', '0', '0', '0']

        # n_hd_cat = 2 (control, hd)
        num_perm = 100

        struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']

        # lets do stats structure by structure
        for k in range(0, len(struct)):

                print "/======================================================="
                print "| Working on the " + struct[k]
                print "\======================================================="
                
                cur_io_dir = io_path + struct[k] + '/'
                # Create the output directory if needed
                if not (os.path.isdir(cur_io_dir)):
                        os.makedirs(cur_io_dir)

                # Create the lpts_file
                lpts = '%slpts_file_%s.txt' %(cur_io_dir, struct[k])
                design = '%sdesign_file.txt' %(cur_io_dir)
                tp = '%stp_file.txt' %(cur_io_dir)
                fixed = '%sfixed_effects_%s.txt' %(cur_io_dir, struct[k])
                random = '%srandom_effects_%s.txt' %(cur_io_dir, struct[k])
                template = '%s%s_time_00.lpts' %(template_file_dir, struct[k])
                out_file = '%shypothesis_ctrl_high.txt' %(cur_io_dir)
                
                # Create test path
                test_path = hypothesis_test_path + ' ' + lpts + ' ' + design + ' ' + tp + ' ' + template + ' ' + fixed + ' ' + random + ' ' + str(num_perm) + ' ' + out_file + ' &'
                os.system(test_path)
                
if __name__ == "__main__":
	main()


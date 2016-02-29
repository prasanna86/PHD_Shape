import os
import math

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	git_path = "/home/sci/prasanna/work/git/PHD_ShapeGrant/deformetrica/tutorial"
        lpts_path = "%s/shape_lpts_for_stats/" %(git_path)
	output_base_path = "%s/mixed_effects_results/" %(git_path)
        mixed_effects_test_path = "/home/sci/prasanna/work/git/PHD_ShapeGrant/shape-lme/bin/mixed-effects-covariate-estimation"

	if not (os.path.isdir(output_base_path)):
		os.system('mkdir ' + output_base_path)

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	times = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]
        sex = ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
        hd = ['3', '3', '3', '0', '0', '3', '0', '3', '3', '3', '0', '0', '0', '0']

	segs = ['left_caudate_time_', 'right_caudate_time_', 'left_putamen_time_', 'right_putamen_time_']
        struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']

        # lets do stats structure by structure
        for k in range(0, len(struct)):
                
                cur_output_dir = output_base_path + struct[k] + "/"
                # Create the output directory if needed
                if not (os.path.isdir(cur_output_dir)):
                        os.makedirs(cur_output_dir)

                # Create the lpts_file
                lpts_file_path = "%slpts_file_%s.txt" %(cur_output_dir, struct[k])
                design_file_path = "%sdesign_file.txt" %(cur_output_dir)
                tp_file_path = "%stp_file.txt" %(cur_output_dir)
		
                lpts_file = open(lpts_file_path,'w')
                design_file = open(design_file_path,'w')
                tp_file = open(tp_file_path,'w')

                for i in range(0, len(subjects)): 
		
                        print "/======================================================="
                        print "| Working on subject " + subjects[i]
                        print "\======================================================="

                        num_cur_timepts = len(times[i])
                        min_time = float(times[i][0])
                        time_range = float(times[i][num_cur_timepts-1]) - float(times[i][0])
                        timepts_to_copy = []
	
                        # We know we want the first and last time points, 
                        # but also intermediate, which we must calculate
                        for j in range(0, num_cur_timepts):

                                time_pt = int(round(((float(times[i][j]) - min_time) / time_range)*29.0))
                                timepts_to_copy.append(time_pt)

                        cur_lpts_dir = lpts_path + subjects[i] + "/"

                        for j in range(0, len(timepts_to_copy)):

                                # Copy the data we need over from the lpts directory
                                subj_timept_lpts_file = "%s%s%0.2d.lpts\n" %(cur_lpts_dir, segs[k], j)
                                lpts_file.write(subj_timept_lpts_file)

                                if(hd[i] == str(0)):
                                        design_matrix_row = str(sex[i]) + " 0 0 0 " + str(times[i][j]) + "\n"
                                elif(hd[i] == str(1)):
                                        design_matrix_row = str(sex[i]) + " 1 0 0 " + str(times[i][j]) + "\n"
                                elif(hd[i] == str(2)):
                                        design_matrix_row = str(sex[i]) + " 0 1 0 " + str(times[i][j]) + "\n"
                                else:
                                        design_matrix_row = str(sex[i]) + " 0 0 1 " + str(times[i][j]) + "\n"

                                design_file.write(design_matrix_row)
                                
                        tp_row = str(len(timepts_to_copy)) + "\n"
                        tp_file.write(tp_row)
                
                lpts_file.close()
                design_file.close()
                tp_file.close()

                mixed_effects_test_call_string = mixed_effects_test_path + ' ' + lpts_file + ' ' + design_file + ' ' + tp_file + ' ' + str(0)

if __name__ == "__main__":
	main()


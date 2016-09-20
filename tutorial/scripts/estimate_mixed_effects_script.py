import os
import math

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	input_path = '../shape_lpts_for_stats/'
        output_path = '../mixed_effects_results/'
        me_test_path = '../shape-lme/bin/mixed-effects-covariate-estimation'

	if not (os.path.isdir(output_path)):
		os.system('mkdir ' + output_path)

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	times = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]
        sex = ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
        hd = ['3', '3', '3', '0', '0', '3', '0', '3', '3', '3', '0', '0', '0', '0']

        n_sex = 1
        n_hd_cat = 2

        segs = ['left_caudate_time_', 'right_caudate_time_', 'left_putamen_time_', 'right_putamen_time_']
        struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']

        # lets do stats structure by structure
        for k in range(0, len(struct)):

                print "/======================================================="
                print "| Working on the " + struct[k]
                print "\======================================================="
                
                cur_output_dir = output_path + struct[k] + '/'
                # Create the output directory if needed
                if not (os.path.isdir(cur_output_dir)):
                        os.makedirs(cur_output_dir)

                # Create the lpts_file
                lpts_file_path = '%slpts_file_%s.txt' %(cur_output_dir, struct[k])
                design_file_path = '%sdesign_file.txt' %(cur_output_dir)
                tp_file_path = '%stp_file.txt' %(cur_output_dir)
                out_fixed_effects = '%sfixed_effects_%s.txt' %(cur_output_dir, struct[k])
                out_random_effects = '%srandom_effects_%s.txt' %(cur_output_dir, struct[k])
                template_file = '%s%s/%s_00.lpts' %(input_path, subjects[3], segs[k])

                lpts_file = open(lpts_file_path,'w')
                design_file = open(design_file_path,'w')
                tp_file = open(tp_file_path,'w')
                
                for i in range(0, len(subjects)): 
		
                        print "/======================================================="
                        print "         | Working on subject " + subjects[i]
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

                        for j in range(0, len(timepts_to_copy)):

                                # Copy the data we need over from the lpts directory
                                subj_timept_lpts_file = '%s%s/%s%0.2d.lpts\n' %(input_path, subjects[i], segs[k], j)
                                lpts_file.write(subj_timept_lpts_file)

                                ia_tag = str(1) + " " + str(times[i][j])

                                if(n_sex == str(2)):
                                        sex_tag = str(" ") + str(sex[i]) + " " + str(times[i][j] * sex[i]) + str(" ")
                                else:
                                        sex_tag = str(" ")
                                 
                                if(n_hd_cat == 2 and hd[i] == str(0)):
                                        hd_tag = str(0) + " " + str(0)
                                elif(n_hd_cat == 2 and hd[i] != str(0)):
                                        hd_tag = str(1) + " " + str(times[i][j])
                                elif(n_hd_cat == 3 and hd[i] == str(0)):
                                        hd_tag = str(0) + " " + str(0) + " " + str(0) + " " + str(0)
                                elif(n_hd_cat == 3 and hd[i] == str(1)):
                                        hd_tag = str(1) + " " + str(times[i][j]) + " " + str(0) + " " + str(0)
                                elif(n_hd_cat == 3 and hd[i] == str(2)):
                                        hd_tag = str(0) + " " + str(0) + " " + str(1) + " " + str(times[i][j])
                                elif(n_hd_cat == 4 and hd[i] == str(0)):
                                        hd_tag = str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0)
                                elif(n_hd_cat == 4 and hd[i] == str(1)):
                                        hd_tag = str(1) + " " + str(times[i][j]) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0)
                                elif(n_hd_cat == 4 and hd[i] == str(2)):
                                        hd_tag = str(0) + " " + str(0) + " " + str(1) + " " + str(times[i][j]) + " " + str(0) + " " + str(0)
                                elif(n_hd_cat == 4 and hd[i] == str(3)):
                                        hd_tag = str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(1) + " " + str(times[i][j])
                                else:
                                        hd_tag = ""
                                d_row = ia_tag + sex_tag + hd_tag + "\n"
                                design_file.write(d_row)
                                
                        tp_row = str(len(timepts_to_copy)) + "\n"
                        tp_file.write(tp_row)
                
                lpts_file.close()
                design_file.close()
                tp_file.close()

                os.system(me_test_path + ' ' + lpts_file_path + ' ' + design_file_path + ' ' + tp_file_path + ' ' + template_file + ' ' + out_fixed_effects + ' ' + out_random_effects)
                
                
if __name__ == "__main__":
	main()


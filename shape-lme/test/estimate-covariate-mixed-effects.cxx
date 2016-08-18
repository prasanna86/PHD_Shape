#include "primitives.hxx"
#include <fstream>
#include <iostream>

MatrixType procrustes(MatrixType & p, MatrixType & q)
{
  MatrixType m = p * strans(q);

  MatrixType U, V; VectorType s;
  svd_econ(U, s, V, m, 'b');

  MatrixType rotation = U * strans(V);
  MatrixType qRot = rotation * q;
  return qRot;
}

void runTests(MatrixType & X, MatrixType & design, TimePointsType & tp, MatrixType & fixed, MatrixType & random)
{
  // cout << "Estimating params using ME" << endl;
  VectorType expl; expl.zeros(design.n_rows);
  expl = design.col(1);

  unsigned int num_fixed = design.n_cols;

  // Number of samples
  unsigned int nsh = static_cast<double>(X.n_cols);
  unsigned int nr = X.n_rows; //number of points*3
  unsigned int nsub = tp.size();

  MatrixType identity_2;
  identity_2.eye(2,2);
  
  MatrixType * Ws = NULL, * Vs = NULL, * identity_n = NULL, * Xp = NULL, * Zp = NULL;
  Ws = new MatrixType[nsub];
  Vs = new MatrixType[nsub];
  identity_n = new MatrixType[nsub];
  Xp = new MatrixType[nsub];  // Xp: Design matrix. for fixed part
  Zp = new MatrixType[nsub];  // Zp: Design matrix. for random part
  
  VectorType * y = NULL, * residual = NULL;
  y = new VectorType[nsub];
  residual = new VectorType[nsub];

  int obs_thus_far = 0;

  for(int k = 0; k < nsub; k++)
  {
    Vs[k].set_size(tp(k), tp(k));
    Ws[k].set_size(tp(k), tp(k));
    
    identity_n[k].eye(tp(k), tp(k));
    Xp[k].set_size(tp(k), num_fixed);
    Zp[k].set_size(tp(k), 2);

    y[k].zeros(tp(k));
    residual[k].zeros(tp(k));
  
    for (int l = 0; l < tp(k); l++)
    {
      Zp[k](l, 0) = 1;
      Zp[k](l, 1) = expl(obs_thus_far + l);

      Xp[k](l, 0) = 1;
      for(int m = 1; m < num_fixed; m++)
      {
	Xp[k](l, m) = design(obs_thus_far + l, m);
      }
    }
 
    obs_thus_far += tp(k);
  }
  
  MatrixType sum_mat1; sum_mat1.zeros(num_fixed, num_fixed);
  VectorType sum_mat2; sum_mat2.zeros(num_fixed);
  double ecorr = 0.0;
  double tracevar = 0.0;
  MatrixType bscorr; bscorr.zeros(2, 2);
  MatrixType bsvar; bsvar.zeros(2, 2);
  VectorType tempvect_fixed; tempvect_fixed.zeros(num_fixed);
  VectorType tempvect_rand; tempvect_rand.zeros(num_fixed);

  double sigma2s = 0, sigma2s_new = 1; // variance of error
  MatrixType Ds, Ds_new; // covariance matrix of random parameters (2x2)
  Ds.zeros(2,2); Ds_new.eye(2,2);
  double epsilon = 1.0e-05;
  double sigmaDiff = 2 * epsilon, Ds_norm_diff = 2 * epsilon;
  unsigned int em = 0;

  for(int i = 0; i < nr; i++) //for all points (x,y,z coordinates)
  {
    em = 0;
    obs_thus_far = 0;

    for(int k = 0; k < nsub; k++)
    {
      y[k].fill(0.0);
  
      for (int l = 0; l < tp(k); l++)
      	y[k](l) = X(i, obs_thus_far + l);
 
      obs_thus_far += tp(k);
    }

    while(em < 500 && (sigmaDiff > epsilon || Ds_norm_diff > epsilon)) // EM iterations
    {
      sigma2s = sigma2s_new;
      Ds = Ds_new;

      sum_mat1.fill(0.0); sum_mat2.fill(0.0);
      ecorr = 0.0;
      tracevar = 0.0;
      bscorr.fill(0.0);
      bsvar.fill(0.0);
     
      for (int k = 0; k < nsub; k++)
      {
	Vs[k] = (identity_n[k] * sigma2s) + Zp[k] * Ds * strans(Zp[k]);
        Ws[k] = inv(Vs[k]);
       
        sum_mat1 += strans(Xp[k]) * Ws[k] * Xp[k];
        sum_mat2 += strans(Xp[k]) * Ws[k] * y[k];
      }

      tempvect_fixed = solve(sum_mat1, sum_mat2);
      fixed.col(i) = tempvect_fixed;

      for (int k = 0; k < nsub; k++)
      {
        tempvect_rand = Ds * strans(Zp[k]) * Ws[k] * (y[k] - (Xp[k] * fixed.col(i)));
        random.col(i * nsub + k) = tempvect_rand;
	residual[k] = y[k] - (Xp[k] * fixed.col(i)) - (Zp[k] * random.col(i * nsub + k));	
        ecorr += dot(residual[k], residual[k]);
        tracevar += (tp(k) - sigma2s * trace(Ws[k]));
        bscorr += random.col(i * nsub + k) * strans(random.col(i * nsub + k));
        bsvar += (identity_2 - (strans(Zp[k]) * Ws[k] * Zp[k] * Ds));
      }
      
      sigma2s_new = (ecorr + sigma2s * tracevar) / static_cast<double>(nsh);
      Ds_new = (bscorr + Ds * bsvar) / static_cast<double>(nsub);

      sigmaDiff = abs(sigma2s_new - sigma2s);
      Ds_norm_diff = norm(Ds_new - Ds, "fro");

      em += 1;
    } //endfor EM iterations

    sigmaDiff = 2 * epsilon;
    Ds_norm_diff = 2 * epsilon;

    // cout << "particle, em iter, Ds_new norm, sigma2 = " << i+1 << " " << em << " " << norm(Ds_new, "fro") << " " << sigma2s_new << endl;
  } //endfor all points on shape (x,y & z)
  
  delete [] Vs;
  delete [] Ws;
  delete [] identity_n;
  delete [] Xp;
  delete [] Zp;
  delete [] y;
  delete [] residual;
}

int main(int argc, char ** argv)
{
  // if(argc != 5)
  // {
  //   cout << "Usage: test/estimate-covariate-mixed-effects param_file template_file out_fixed out_random"
  // 	 << endl;
  //   return EXIT_FAILURE;
  // }
  
  // variables needed to read param file
  ifstream input_file, shape_file;
  ScalarType xVal, px, py, pz;
  string filename;
  unsigned int gid; bool s;

  // Variables we will need for regression
  unsigned int subjectNumber, prevSubjectNumber;

  // The Lists
  std::vector<unsigned int> subj_nums;
  std::vector<unsigned int> groupIDs;
  std::vector<bool> sex_of_subjs;
  ScalarListType x;

  input_file.open(argv[1]);
  input_file >> subjectNumber >> gid >> s >> xVal >> filename;
  prevSubjectNumber = subjectNumber;

  shape_file.open(filename.c_str());
  cout << "Reading file: " << " " << filename << endl;
  std::string line;
  while (std::getline(shape_file, line))
    ++np;
  shape_file.close();
  
  unsigned int index = 0, t = 0, j = 0;
  unsigned int np, nsh, nsub, num_fixed, num_cat;

  // reading param file
   while(!input_file.eof())
   {
     while(subjectNumber == prevSubjectNumber)
     {
       // incrementing time point
       t++;

       groupIDs.push_back(gid);
       sex_of_subjs.push_back(s);
       x(j) = xVal;
       subj_nums.push_back(subjectNumber);
       j += 1;

       inFile >> subjectNumber >> gid >> s >> xVal >> filename;
       if(inFile.fail())
	 break;
     }
     tp[index] = t; ++index;
     // re-intializing
     t = 0;
     prevSubjectNumber = subjectNumber;
   }

   MatrixType X, Y; X.zeros(3*np, nsh); Y.zeros(3*np, nsh);
   MatrixType design; design.zeros(nsh, num_fixed);
   TimePointsType tp; tp.set_size(nsub);



  // // loading design file
  // unsigned int num_fixed = 0;
  // design_file.open(argv[2]);
  // std::getline(design_file, line);
  // std::istringstream iss(line);
  // do
  // {
  //   std::string sub;
  //   iss >> sub;
  //   if (sub.length())
  //     ++num_fixed;
  // }
  // while(iss);
  // design_file.close();

  // // std::cout << "num_fixed: " << num_fixed << "\n";

  // // counting number of subjects from timepts file
  // unsigned int nsub = 0;
  // timepoints_file.open(argv[3]);
  // while (std::getline(timepoints_file, line))
  //   ++nsub;

  // timepoints_file.close();
  // // std::cout << "nsub: " << nsub << "\n";

  // // counting number of particles
  // unsigned int np = 0;
  // lpts_paramfile.open(argv[1]);
  // lpts_paramfile >> filename;
  // shapefile.open(filename.c_str());
  // while (std::getline(shapefile, line))
  //   ++np; 
  // // std::cout << "np: " << np << "\n";
  // shapefile.close();
  // lpts_paramfile.close();

  // MatrixType X, Y; X.zeros(3*np, nsh); Y.zeros(3*np, nsh);
  // MatrixType design; design.zeros(nsh, num_fixed);
  // TimePointsType tp; tp.set_size(nsub);

  // lpts_paramfile.open(argv[1]);

  // unsigned int j = 0, k = 0;
  // while(lpts_paramfile >> filename)
  // {
  //   shapefile.open(filename.c_str());
  //   // cout << "Reading file: " << filename << endl;
  //   int i = 0;
  //   while(shapefile >> px >> py >> pz)
  //   {
  //     X(3 * i, j) = px;
  //     X(3 * i + 1, j) = py;
  //     X(3 * i + 2, j) = pz;
  //     i += 1;
  //   }
  //   shapefile.close();
  //   ++j;
  // }

  // lpts_paramfile.close();

  // // loading design matrix
  // j = 0;
  // design_file.open(argv[2]);
  // while(design_file >> entry)
  // {
  //   design(j, k) = entry;
  //   ++k;
  //   if(k == num_fixed)
  //   {
  //     ++j;
  //     k = 0;
  //   }
  // }

  // design_file.close();

  // // loading timepoints
  // j = 0;
  // timepoints_file.open(argv[3]);
  // while(timepoints_file >> t)
  // {    
  //   tp[j] = t;
  //   j += 1;
  // }

  // timepoints_file.close();

  // // parameters of interest
  // MatrixType fixed; // slopes + intercepts for all points
  // MatrixType random; // slopes + intercepts for all groups, for all points
  // fixed.set_size(num_fixed, 3 * np);
  // random.set_size(2, 3 * np * nsub);

  // // align shapes to template
  // MatrixType template_shape; template_shape.zeros(3, np);
  // // load in a nice template shape
  // shapefile.open(argv[4]);
  // j = 0;
  // while(shapefile >> px >> py >> pz)
  // {
  //   template_shape(0, j) = px;
  //   template_shape(1, j) = py;
  //   template_shape(2, j) = pz;
  //   j += 1;
  // }
  // shapefile.close();

  // // cout << "template = " << endl << strans(template_shape) << endl;

  // // removing translations from template
  // VectorType mean_xyz; mean_xyz.zeros(3); 
  // for(int k = 0; k < np; k++)
  //   mean_xyz += template_shape.col(k);

  // mean_xyz /= static_cast<double>(np);

  // for(int k = 0; k < np; k++)
  //   template_shape.col(k) = template_shape.col(k) - mean_xyz;

  // j = 0; // reinitializing

  // // ...and align everything to the template
  // for(int k = 0; k < nsh; k++)
  // {
  //   MatrixType current_shape = X.col(k);
  //   current_shape.reshape(3, np);

  //   // removing translations from current shape
  //   mean_xyz.zeros(3); 
  //   for(int k = 0; k < np; k++)
  //     mean_xyz += current_shape.col(k);

  //   mean_xyz /= static_cast<double>(np);

  //   for(int k = 0; k < np; k++)
  //     current_shape.col(k) = current_shape.col(k) - mean_xyz;

  //   // procrustes align to template
  //   current_shape = procrustes(template_shape, current_shape);
    
  //   // reshape back
  //   current_shape.reshape(3 * np, 1);
  //   Y.col(k) = current_shape;
  // }
    
  // // estimate model with aligned shapes
  // runTests(Y, design, tp, fixed, random);

  // // output fixed params
  // ofstream out;
  // out.open(argv[5]);
 
  // for(int i = 0; i < np * 3; i++)
  // {
  //   for(int j = 0; j < num_fixed; j++)
  //     out << fixed(j, i) << " "; 
  //   out << endl;
  // }

  // out.close();
  // // output random params
  // out.open(argv[6]);

  // for(int i = 0; i < nsub * np * 3; i++)
  // {
  //   for(int j = 0; j < 2; j++)
  //     out << random(j, i) << " ";
  //   out << endl;
  // }

  // out.close();
  return EXIT_SUCCESS;
}

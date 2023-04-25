void write_data (Params &par, Operators &opr, int i)
{
    // Writing data into a file in the format of:
    // index, density, real potential.
  double density[(int)pow(opr.size,2)];
  
  for (size_t j = 0; j < pow(opr.size,2); ++j) 
  {
    density[j] = pow(abs(opr.wfc[j]), 2); // norm(wfc)**2
  }

  if (par.im_time)
  {
      double sum = 0;

      for (size_t j = 0; j < pow(opr.size,2); ++j)
      {
	  sum += density[j];
      }

      sum *= par.dx;

      for (size_t j = 0; j < pow(opr.size,2); ++j)
      {
	  opr.wfc[j] /= sqrt(sum);
      }
  }
    
  // file path
  std::stringstream filename_fdm, filename_particle ;
  filename_fdm << "../output/fdm/output" << i << ".dat";
  filename_particle << "../output/particle/output" << i << ".dat";

  // file name
  std::ofstream ffdm = std::ofstream(filename_fdm.str());
  std::ofstream fparticle = std::ofstream(filename_particle.str());

  if (ffdm)
  {
      for (int i = 0; i < opr.size; ++i)
      {
	for (int j = 0; j < opr.size; ++j)
	{
	  std::stringstream data_fdm;
	  data_fdm << par.x[i] << "\t" << par.y[j] << "\t" << density[j+i*opr.size] << "\t" << opr.phi[j+i*opr.size] << "\t" << opr.ax[j+i*opr.size] << "\t" << opr.ay[j+i*opr.size] <<"\n";
	  ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
//	  if (i == 0) {std::cout << opr.wfc[j] << "\t";}
	}
      }
  }
  ffdm.close();

  if (fparticle)
  {
    for (int j = 0; j < par.np*2; ++j)
    {
      std::stringstream data_particle;
      data_particle << par.px[j] << "\t" << par.py[j] << "\t" <<
  	 par.pvx[j] <<"\t" << par.pvy[j] <<"\n";
      fparticle.write(data_particle.str().c_str(), data_particle.str().length());
    }
  }
  fparticle.close();
}

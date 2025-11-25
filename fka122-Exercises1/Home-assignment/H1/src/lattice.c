void init_fcc(double **positions, int N, double lattice_param)
{
    double base[4][3] = {
	{0., 0., 0.},
	{0.5, 0.5, 0.},
	{0., 0.5, 0.5},
	{0.5, 0., 0.5}
    };
    int ind = 0;
    for (int x = 0; x < N; x++){
        for (int y = 0; y < N; y++){
            for (int z = 0; z < N; z++){
		for(int b = 0; b < 4; b++){
		    positions[ind][0] = lattice_param * (base[b][0] + x);
                    positions[ind][1] = lattice_param * (base[b][1] + y);
                    positions[ind][2] = lattice_param * (base[b][2] + z);
		    ind += 1;
		}
            }
        }
    }
}

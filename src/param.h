struct param // parameter struct to be passed to GSL ODE integrators
{
    int n_iso;  // number of isotopes included in network
    double T;   // temperature
    double rho; // mass density
};

//  Darcy flow solver

{
    int counter=0;
    while(counter < nCycle)
    {
        counter++;
        double residual = solve
        ( fvm::laplacian(H*H*H,p) ).initialResidual();

        if (DEBUG)
        Info << "Iteration "    << counter 
             << ": residual = " << residual << endl;

        if (residual < tolerance)
        {
            Info << "Darcy solver: Converged in "
                 << counter << " steps" << endl;
            break;                                  // Done
        }                
        else if (counter == nCycle)
        {
            Info << "Darcy: FAILED to converge in "
                 << nCycle << " steps" << nl << endl;
            exit(1);                                // No convergence 
        }
    }
}

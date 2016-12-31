//  Nonlinear reaction kinetics + transport correction  //

    double stoiTable[1][2] = {{1.0,-1.0}};      // Stoichiometric coeffs

    {
        int iter, maxIter=100;
        double Sh=8.0, p=4.0, Ct=0.2, len=0.02; // Parameters
        double c, cR, x, g;
        double w, w1, w2, r, dr;
        volScalarField& C = CList[0]();

        R  = 2*k0;
        dR = 2*k0;
        forAll (C, i)
        {
            g  = 2*(k0/D0).value()*H[i]/Sh;     // Transport correction
            c  = C[i]/(1+g);                    // Linear initial guess
            if (c < Ct)  c = Ct;                // Ct as initial min val
            for (iter = 0; iter < maxIter; iter++)
            {
                cR = c;                         // Current c as refValue
                w  = std::tanh((cR-Ct)/len);
                w1 = 0.5*(1+w);
                w2 = 0.5*(1-w);
                x  = std::pow(cR/Ct,p-1);
                r  = (w1 + w2*x)*cR;
                dr = (w1 + w2*x*p) + (1-w*w)*(1-x)*cR/(2*len);
                c  = cR - (cR+g*r-C[i])/(1+g*dr);
                if (fabs(c-cR) < tolerance)  break;
            }
            if (iter == maxIter)
            {
                Info << nl
                     << "reacModel: FAILED to converge in "
                     << maxIter << " iterations" << nl << endl;
                     exit(1);                   // No convergence
            }
            R[i]  *= r;
            dR[i] *= dr;
        }
         R.correctBoundaryConditions();         // Update proc patches
        dR.correctBoundaryConditions();

        RList[0]()  == R;
        dRList[0]() == dR;
    }
